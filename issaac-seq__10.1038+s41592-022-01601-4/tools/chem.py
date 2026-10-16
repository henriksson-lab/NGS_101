"""Source-backed ISSAAC-seq plate-workflow chemistry."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import (Construct, MolecularState, Row, Scene, Workflow,
                      duplex_rows, feature)
from multimodal_spatial import finish_pcr_branch
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind
from tagmentation import loaded_tn5_adapter

TITLE = "ISSAAC-seq — accessibility and RNA from the same cell"
NOTES = "01_issaac-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/s41592-022-01601-4">'
          'Xu et al., <i>Nature Methods</i> (2022)</a>.')
SUMMARY = ("One S5/S7 transposome tags accessible chromatin, reverse transcription makes "
           "RNA/DNA hybrids, and a second transposome tags those hybrids: S7/S7 in the "
           "plate workflow drawn here, S5/S5 in the droplet workflow. Distinct primer "
           "pairs then recover ATAC and RNA libraries from each cell.")
CAVEAT = ("The exact constructs shown are the plate workflow from Supplementary Table 4. "
          "Published Supplementary Protocols 1 and 2 specify S5/S5 for the second "
          "droplet tagmentation and S7/S7 for the second plate tagmentation. The "
          "droplet workflow also replaces the plate indexes with a 16-nt 10x scATAC "
          "gel-bead cell barcode and uses the published TruseqR2 capture/P7 RNA branch.")

ME_S5 = nx.ADAPTOR_S5
ME_S7 = nx.ADAPTOR_S7
ME_BOTTOM = nx.ME_RC
TRUSEQ_R1_SHORT = "CTACACGACGCTCTTCCGATCT"
TRUSEQ_R2_SHORT = "CAGACGTGTGCTCTTCCGATCT"
UMI_LENGTH = 10
DT_LENGTH = 32
PLATE_READ_LENGTHS = {"Read 1": 150, "Index 1 (i7)": 8,
                      "Index 2 (i5)": 8, "Read 2": 150}
ATAC_P5_PRIMER = il.P5 + "N" * 8 + nx.S5
RNA_P5_PRIMER = il.P5 + "N" * 8 + TRUSEQ_R1_SHORT
P7_PRIMER = il.P7 + "N" * 8 + nx.S7


def tn5_adapter_scene(entry: str) -> Scene:
    if entry == "S5":
        top = [seg("S5", nx.S5, "s5"), seg("mosaic end", nx.ME, "me")]
    elif entry == "S7":
        top = [seg("S7", nx.S7, "s7"), seg("mosaic end", nx.ME, "me")]
    else:
        raise ValueError("ISSAAC Tn5 adapter is S5 or S7")
    return loaded_tn5_adapter(
        top, mosaic_end="mosaic end", label=f"ME_{entry}",
        bottom_name="ME_bottom",
        bottom_label="5′-phosphorylated, phosphorothioated ME_bottom",
        bottom_sequence=ME_BOTTOM, mod5="p", mod3="inverted dT")


def atac_tagged() -> Construct:
    return Construct([
        seg("S5", nx.S5, "s5"), seg("left mosaic end", nx.ME, "me"),
        seg("accessible chromatin", "X" * 36, placeholder=True),
        seg("right mosaic end reverse complement", nx.ME_RC, "me"),
        seg("S7 reverse complement", nx.S7_RC, "s7"),
    ], name="S5/S7-tagged accessible chromatin")


ATAC_TAGGED = atac_tagged()


def rt_scene() -> Scene:
    mrna = [seg("transcript", "X" * 34, placeholder=True),
            seg("poly(A)", "A" * DT_LENGTH)]
    primer = [seg("TruSeq R1 short arm", TRUSEQ_R1_SHORT, "r1"),
              seg("UMI", "U" * UMI_LENGTH, "umi", placeholder=True,
                  feature=feature("umi", "umi", "random")),
              seg("oligo(dT)32", "T" * DT_LENGTH),
              seg("VN anchor", "VN", placeholder=True)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="polyadenylated RNA")
    sc.anneal("capture", primer, to="mRNA", pair=("oligo(dT)32", "poly(A)"),
              label="5′-biotin TruseqR1_oligo_dT",
              unpaired=("VN anchor",), mod5="biotin")
    sc.arrow("capture", "reverse transcription")
    sc.mark("capture", "UMI", "10-nt molecule identity")
    sc.labels("capture")
    return sc


def rna_tagged() -> Construct:
    return Construct([
        seg("TruSeq R1 short arm", TRUSEQ_R1_SHORT, "r1"),
        seg("UMI", "U" * UMI_LENGTH, "umi", placeholder=True,
            feature=feature("umi", "umi", "random")),
        seg("copied oligo(dT)32", "T" * DT_LENGTH),
        seg("VN anchor", "VN", placeholder=True),
        seg("3-prime cDNA", "X" * 36, placeholder=True),
        seg("mosaic end reverse complement", nx.ME_RC, "me"),
        seg("S7 reverse complement", nx.S7_RC, "s7"),
    ], name="S7-tagged 3-prime RT product")


RNA_TAGGED = rna_tagged()


def _indexes():
    return (
        seg("i5 plate/modality index", "J" * 8, "cbc", placeholder=True,
            feature=feature("plate_modality_index", "sample_index", "fixed")),
        seg("i7 well barcode reverse complement", "I" * 8, "cbc", placeholder=True,
            feature=feature("well_cell_barcode", "cell_barcode", "fixed")),
    )


def final_atac() -> Construct:
    i5, i7 = _indexes()
    return Construct([seg("P5", il.P5, "p5"), i5, *list(ATAC_TAGGED), i7,
                      seg("P7 reverse complement", il.P7_RC, "p7")],
                     name="plate ISSAAC-seq ATAC library")


def final_rna() -> Construct:
    i5, i7 = _indexes()
    # The P5 primer extends the short capture arm to the complete TruSeq Read 1 site.
    body = list(RNA_TAGGED)[1:]
    return Construct([seg("P5", il.P5, "p5"), i5,
                      seg("TruSeq Read 1 arm", il.TRUSEQ_READ1, "r1"),
                      *body, i7, seg("P7 reverse complement", il.P7_RC, "p7")],
                     name="plate ISSAAC-seq RNA library")


ATAC = final_atac()
RNA = final_rna()
ATAC_PRIMERS = tuple(sp.NEXTERA[k] for k in ("R1", "I1", "I2", "R2"))
RNA_PRIMERS = (sp.TRUSEQ["R1"], sp.NEXTERA["I1"],
               sp.TRUSEQ["I2"], sp.NEXTERA["R2"])
for _name, _lib, _primers in (("ATAC", ATAC, ATAC_PRIMERS),
                               ("RNA", RNA, RNA_PRIMERS)):
    _problems = sp.verify(_lib, _primers)
    if _problems:
        raise ValueError(f"ISSAAC-seq {_name}: " + "; ".join(_problems))

FINAL_LIBRARIES = (
    ("Plate ATAC library", ATAC, ATAC_PRIMERS,
     "S5/S7-tagged accessible chromatin; i7 is the sorted-well cell identity and i5 identifies the plate/modality.",
     "All four sequencing primers are the standard Nextera set and land by exact sequence."),
    ("Plate RNA library", RNA, RNA_PRIMERS,
     "The capture-derived TruSeq end retains the 10-nt UMI; the S7-tagmented end is Nextera. i7 identifies the sorted well.",
     "This is a hybrid library: TruSeq Read 1 and Index 2 primers bind one end, while Nextera Index 1 and Read 2 primers bind the other."),
)
READ_LENGTHS = {heading: PLATE_READ_LENGTHS for heading, *_ in FINAL_LIBRARIES}
READ_USES = {
    "Plate ATAC library": {"Read 1": "include", "Index 1 (i7)": "extract",
                           "Index 2 (i5)": "demultiplex", "Read 2": "include"},
    "Plate RNA library": {"Read 1": "extract_then_include", "Index 1 (i7)": "extract",
                          "Index 2 (i5)": "demultiplex", "Read 2": "include"},
}


def oligos():
    return (
        f"<b>ME_S5 / ME_S7</b> — 5′-{ME_S5}-3′ / 5′-{ME_S7}-3′",
        f"<b>ME_bottom</b> — 5′-Phos-C*T*G*T*C*T*C*T*T*A*T*A*C*A*C*A*T*C*/iInvdT/-3′",
        f"<b>TruseqR1_oligo_dT</b> — 5′-biotin-{TRUSEQ_R1_SHORT}[10-nt UMI]{'T' * DT_LENGTH}VN-3′",
        f"<b>TruseqR2_oligo_dT (droplet)</b> — 5′-biotin-{TRUSEQ_R2_SHORT}[10-nt UMI]{'T' * DT_LENGTH}VN-3′",
        f"<b>Truseq_S5_short / Nextra_S5_short</b> — 5′-{TRUSEQ_R1_SHORT}-3′ / 5′-{nx.S5}-3′",
        "<b>N7001, example i7 well-index primer</b> — 5′-CAAGCAGAAGACGGCATACGAGATAGTAAACCGTCTCGTGGGCTCGG-3′",
        "<b>RNA_plate_S502, example i5 RNA primer</b> — 5′-AATGATACGGCGACCACCGAGATCTACACCTCTCTATACACTCTTTCCCTACACGACGCTCTTCCG-3′",
        "<b>ATAC_plate_SS502, example i5 ATAC primer</b> — 5′-AATGATACGGCGACCACCGAGATCTACACCTCTCTATTCGTCGGCAGCGTC-3′",
    )


def workflow() -> Workflow:
    initial = Scene.duplex([seg("native chromatin", "X" * 42, placeholder=True)],
                           label="fixed nucleus")
    wf = Workflow(MolecularState("Fixed nuclei containing chromatin and poly(A) RNA",
                                 tuple(initial.rows())))
    wf.react("Anneal ME_S5 and ME_S7 to modified ME_bottom and load S5/S7 Tn5",
             [*tn5_adapter_scene("S5").rows(), *tn5_adapter_scene("S7").rows()],
             name="S5/S7 transposome")
    wf.react("Tagment accessible chromatin in situ",
             duplex_rows(ATAC_TAGGED, label="accessible-DNA branch"),
             name="S5/S7-tagged chromatin")
    wf.react("Prime poly(A) RNA with biotin-TruseqR1_oligo_dT and reverse-transcribe",
             rt_scene().rows(), name="RNA/cDNA hybrid with 10-nt UMI")
    wf.react("Load an S7/S7 transposome and tagment the RNA/cDNA hybrid",
             [*tn5_adapter_scene("S7").rows(),
              *duplex_rows(RNA_TAGGED, label="productive 3′ RNA fragment")],
             name="S7-tagged RT branch",
             note="Only the capture-arm-to-S7 3′ fragment has unlike amplification ends; internal S7/S7 products are suppressed.")
    wf.react("Gap-fill, remove leftover oligos with Exonuclease I, and index-PCR both branches",
             [*duplex_rows(ATAC, label="plate ATAC library"),
              *duplex_rows(RNA, label="plate RNA library")],
             name="Plate ISSAAC-seq libraries",
             note="The i7 well barcode is the cell identity; i5 distinguishes plate and ATAC versus RNA.")
    return wf


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    nuclei = g.input(
        "nuclei", "fixed nuclei containing accessible chromatin and poly(A) RNA",
        MaterialKind.BIOLOGICAL,
        state=Scene.duplex([seg("native chromatin", "X" * 42, placeholder=True)],
                           label="fixed nucleus").rows(), primary=True)

    atac_transposome = g.input(
        "atac_transposome", "S5/S7-loaded chromatin transposome",
        MaterialKind.ADAPTER,
        state=[*tn5_adapter_scene("S5").rows(),
               *tn5_adapter_scene("S7").rows()])
    atac_tn5 = g.input("atac_tn5", "chromatin Tn5 transposase",
                       MaterialKind.CATALYST)
    atac_tagged_nuclei = g.product(
        "atac_tagged_nuclei",
        "nuclei retaining S5/S7-tagged chromatin and unprocessed RNA",
        MaterialKind.MOLECULAR,
        state=[*duplex_rows(ATAC_TAGGED, label="accessible-DNA branch"),
               Row(chunks=[("poly(A) RNA remains available for capture", None, False)])])
    atac_waste = g.product(
        "atac_waste", "non-library chromatin tagmentation products",
        MaterialKind.WASTE)
    g.react(
        "tagment_atac", "Tagment accessible chromatin in situ",
        StepKind.TAGMENTATION,
        ((nuclei, InputRole.TEMPLATE),
         (atac_transposome, InputRole.STRUCTURAL),
         (atac_tn5, InputRole.CATALYST)),
        (atac_tagged_nuclei, atac_waste))
    g.discard(atac_waste,
              "does not retain the selected unlike-ended S5/S7 accessible-DNA fragment")

    capture = g.input(
        "capture_primer", "biotin-TruseqR1-UMI-oligo(dT)32VN primer",
        MaterialKind.OLIGO, state=rt_scene())
    rt = g.input("rt", "reverse transcriptase", MaterialKind.CATALYST)
    copied = g.product(
        "copied", "nuclei retaining tagged chromatin and UMI-labelled RNA/cDNA hybrids",
        MaterialKind.MOLECULAR,
        state=[*duplex_rows(ATAC_TAGGED, label="accessible-DNA branch"),
               *rt_scene().rows()])
    g.react(
        "reverse_transcribe",
        "Prime poly(A) RNA with biotin-TruseqR1-oligo(dT) and reverse-transcribe",
        StepKind.REVERSE_TRANSCRIPTION,
        ((atac_tagged_nuclei, InputRole.TEMPLATE),
         (capture, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)),
        (copied,))

    rna_transposome = g.input(
        "rna_transposome", "S7/S7-loaded RNA/cDNA transposome",
        MaterialKind.ADAPTER, state=tn5_adapter_scene("S7"))
    rna_tn5 = g.input("rna_tn5", "RNA/cDNA-hybrid Tn5 transposase",
                      MaterialKind.CATALYST)
    both_tagged = g.product(
        "both_tagged", "nuclei retaining productive ATAC and RNA fragments",
        MaterialKind.MOLECULAR,
        state=[*duplex_rows(ATAC_TAGGED, label="accessible-DNA branch"),
               *duplex_rows(RNA_TAGGED, label="productive 3-prime RNA fragment")])
    rna_waste = g.product(
        "rna_waste", "suppressed or non-library RNA/cDNA tagmentation products",
        MaterialKind.WASTE)
    g.react(
        "tagment_rna", "Tagment the RNA/cDNA hybrid with S7/S7 Tn5",
        StepKind.TAGMENTATION,
        ((copied, InputRole.TEMPLATE),
         (rna_transposome, InputRole.STRUCTURAL),
         (rna_tn5, InputRole.CATALYST)), (both_tagged, rna_waste))
    g.discard(rna_waste,
              "includes internal S7/S7 products lacking unlike amplification ends")

    fill_polymerase = g.input("fill_polymerase", "gap-fill DNA polymerase",
                              MaterialKind.CATALYST)
    filled = g.product(
        "filled", "gap-filled ATAC and RNA branch templates",
        MaterialKind.MOLECULAR, state=both_tagged.state)
    g.react("gap_fill", "Fill transposition gaps", StepKind.EXTENSION,
            ((both_tagged, InputRole.TEMPLATE),
             (fill_polymerase, InputRole.CATALYST)), (filled,))
    exoi = g.input("exoi", "Exonuclease I", MaterialKind.CATALYST)
    cleaned = g.product(
        "cleaned", "gap-filled templates after free-oligo removal",
        MaterialKind.MOLECULAR, state=filled.state)
    oligo_waste = g.product("oligo_waste", "digested leftover oligos",
                            MaterialKind.WASTE)
    g.react("remove_oligos", "Remove leftover single-stranded oligos with Exonuclease I",
            StepKind.CLEANUP,
            ((filled, InputRole.TEMPLATE), (exoi, InputRole.CATALYST)),
            (cleaned, oligo_waste))
    g.discard(oligo_waste, "is the free-oligo fraction removed before index PCR")

    atac_template = g.product("atac_template", "plate ATAC PCR template",
                              MaterialKind.MOLECULAR, state=ATAC)
    rna_template = g.product("rna_template", "plate RNA PCR template",
                             MaterialKind.MOLECULAR, state=RNA)
    g.react("split", "Allocate ATAC and RNA products to their primer pairs",
            StepKind.SPLIT, ((cleaned, InputRole.TEMPLATE),),
            (atac_template, rna_template), changes_molecule=False)

    reviewed = [
        "nuclei", "atac_transposome", "atac_tn5", "tagment_atac",
        "capture_primer", "rt", "reverse_transcribe", "rna_transposome",
        "rna_tn5", "tagment_rna", "fill_polymerase", "gap_fill", "exoi",
        "remove_oligos", "split",
    ]
    atac_geometry = PCRGeometry.exact(
        ATAC, PCRPrimer("plate ATAC P5/i5/S5", ATAC_P5_PRIMER,
                        anneal_nt=len(nx.S5)),
        PCRPrimer("plate P7/i7/S7", P7_PRIMER, anneal_nt=len(nx.S7)),
        statement="Supplementary Table 4 defines the indexed S5 and S7 plate primers; their 3-prime footprints face inward on the ATAC product.")
    rna_geometry = PCRGeometry.exact(
        RNA, PCRPrimer("plate RNA P5/i5/TruSeq", RNA_P5_PRIMER,
                       anneal_nt=len(TRUSEQ_R1_SHORT)),
        PCRPrimer("plate P7/i7/S7", P7_PRIMER, anneal_nt=len(nx.S7)),
        statement="Supplementary Table 4 defines the hybrid TruSeq/S7 plate-primer pair; their 3-prime footprints face inward on the RNA product.")
    finish_pcr_branch(
        g, "atac", atac_template, ATAC, ATAC_PRIMERS, reviewed,
        label="plate ISSAAC-seq ATAC", read_lengths=PLATE_READ_LENGTHS,
        read_uses=READ_USES["Plate ATAC library"], geometry=atac_geometry)
    finish_pcr_branch(
        g, "rna", rna_template, RNA, RNA_PRIMERS, reviewed,
        label="plate ISSAAC-seq RNA", read_lengths=PLATE_READ_LENGTHS,
        read_uses=READ_USES["Plate RNA library"], geometry=rna_geometry)
    g.declare_reviewed_scope(reviewed, SOURCE)
    return g.require_publishable()
