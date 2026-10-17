"""Source-backed ISSAAC-seq plate-workflow chemistry."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import (Construct, MolecularState, Row, Scene, Segment, Workflow,
                      duplex_rows, feature)
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
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
          "gel-bead cell barcode and uses the published TruseqR2 capture/P7 RNA branch. "
          "Model correction (2026-10-17): Supplementary Protocol 2 adds the i7 well barcode "
          "in-well, by 4 RNA pre-PCR cycles (Truseq_S5_short + N7xxx) and then 6 ATAC+RNA cycles "
          "(+ Nextra_S5_short), before pooling the plate, Exonuclease I digestion and one 8-cycle "
          "final PCR containing ATAC_plate_SS5xx, RNA_plate_S5xx and illumina P7. The earlier page "
          "showed a single split dual-index PCR. The pre-sort gap fill-in uses Maxima H Minus RT "
          "together with Exonuclease I in nuclei (steps 37–39).")

ME_S5 = nx.ADAPTOR_S5
ME_S7 = nx.ADAPTOR_S7
ME_BOTTOM = nx.ME_RC
TRUSEQ_R1_SHORT = "CTACACGACGCTCTTCCGATCT"
TRUSEQ_R2_SHORT = "CAGACGTGTGCTCTTCCGATCT"
UMI_LENGTH = 10
DT_LENGTH = 30  # MOESM4 TruseqR1_oligo_dT / TruseqR2_oligo_dT: T30VN
PLATE_READ_LENGTHS = {"Read 1": 150, "Index 1 (i7)": 8,
                      "Index 2 (i5)": 8, "Read 2": 150}
# Supplementary Table 4 plate primers, with the 8-nt index window as N.
# ATAC_plate_SS5xx = P5 + i5 + Nextra_S5_short; RNA_plate_S5xx = P5 + i5 +
# ACACTCTTTCCCTACACGACGCTCTTCCG (a TruSeq Read 1 prefix ending 4 nt before the
# capture arm's 3' ATCT); N7xxx = Illumina P7 + i7 + S7.  The in-well pre-PCR uses
# Truseq_S5_short / Nextra_S5_short with N7xxx; the final PCR uses the i5 primers
# with the bare Illumina P7 primer (Supplementary Protocol 2 steps 7-9, 48-62).
ATAC_P5_PRIMER = il.P5 + "N" * 8 + nx.S5
RNA_P5_PRIMER = il.P5 + "N" * 8 + il.TRUSEQ_READ1[:-4]
P7_PRIMER = il.P7 + "N" * 8 + nx.S7
ILLUMINA_P7_PRIMER = il.P7
NEXTRA_S5_SHORT = nx.S5


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


def capture_primer() -> list:
    """MOESM4 TruseqR1_oligo_dT, 5'->3' (5'-biotin): handle + N10 + T30 + VN."""
    return [seg("TruSeq R1 short arm", TRUSEQ_R1_SHORT, "r1"),
            seg("UMI", "U" * UMI_LENGTH, "umi", placeholder=True,
                feature=feature("umi", "umi", "random")),
            seg("oligo(dT)30", "T" * DT_LENGTH),
            seg("VN anchor", "VN", placeholder=True)]


def rt_scene() -> Scene:
    mrna = [seg("transcript", "X" * 34, placeholder=True),
            seg("poly(A)", "A" * DT_LENGTH)]
    primer = capture_primer()
    sc = Scene()
    sc.strand("mRNA", mrna, label="polyadenylated RNA")
    sc.anneal("capture", primer, to="mRNA", pair=("oligo(dT)30", "poly(A)"),
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
        seg("copied oligo(dT)30", "T" * DT_LENGTH),
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


def well_preamplicon(tagged: Construct, name: str) -> Construct:
    """In-well pre-PCR product: N7xxx adds the i7 well barcode and P7, not i5."""
    _, i7 = _indexes()
    return Construct([*list(tagged), i7, seg("P7 reverse complement", il.P7_RC, "p7")],
                     name=name)


def _index_window(index: Segment, name: str) -> Segment:
    """An 8-nt index window in primer orientation, keeping the library feature."""
    return seg(name, "N" * 8, "cbc", placeholder=True, feature=index.feature)


def short_primer(name: str) -> Construct:
    """Supplementary Table 4 Truseq_S5_short or Nextra_S5_short primer."""
    if name == "Truseq_S5_short":
        return Construct([seg(name, TRUSEQ_R1_SHORT, "r1")], name=name)
    if name == "Nextra_S5_short":
        return Construct([seg(name, NEXTRA_S5_SHORT, "s5")], name=name)
    raise ValueError(f"unknown ISSAAC short primer: {name!r}")


def n7_primer() -> Construct:
    """Supplementary Table 4 N7xxx family: Illumina P7 + 8-nt i7 well index + S7."""
    _, i7 = _indexes()
    return Construct([seg("P7", il.P7, "p7"), _index_window(i7, "i7 well index"),
                      seg("S7", nx.S7, "s7")], name="N7xxx i7 well-index primer")


def plate_i5_primer(modality: str) -> Construct:
    """Supplementary Table 4 ATAC_plate_SS5xx / RNA_plate_S5xx: P5 + 8-nt i5 + arm."""
    i5, _ = _indexes()
    if modality == "atac":
        arm, name = seg("Nextra S5", nx.S5, "s5"), "ATAC_plate_SS5xx"
    elif modality == "rna":
        arm = seg("TruSeq Read 1 prefix", il.TRUSEQ_READ1[:-4], "r1")
        name = "RNA_plate_S5xx"
    else:
        raise ValueError(f"unknown ISSAAC modality: {modality!r}")
    return Construct([seg("P5", il.P5, "p5"), _index_window(i5, "i5 plate/modality index"),
                      arm], name=name)


ATAC_PRE = well_preamplicon(ATAC_TAGGED, "i7-barcoded ATAC pre-amplicon")
RNA_PRE = well_preamplicon(RNA_TAGGED, "i7-barcoded RNA pre-amplicon")
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
        f"<b>Truseq_S5_short / Nextra_S5_short (in-well pre-PCR)</b> — 5′-{TRUSEQ_R1_SHORT}-3′ / 5′-{NEXTRA_S5_SHORT}-3′",
        f"<b>illumina P7 (final PCR)</b> — 5′-{ILLUMINA_P7_PRIMER}-3′",
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
    wf.react("Gap-fill with Exonuclease I, sort one nucleus per well, pre-PCR with i7, pool, Exonuclease I, then i5 PCR",
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
        state=(tn5_adapter_scene("S5"), tn5_adapter_scene("S7")))
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
        "capture_primer", "biotin-TruseqR1-UMI-oligo(dT)30VN primer",
        MaterialKind.OLIGO,
        state=Construct(capture_primer(), name="5′-biotin TruseqR1_oligo_dT"))
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

    # Supplementary Protocol 2 steps 37-39: one in-nucleus reaction with Maxima H Minus
    # RT, dNTPs and Exonuclease I (37 C, 15 min), before sorting.
    fill_rt = g.input("fill_rt", "Maxima H Minus reverse transcriptase (gap fill-in)",
                      MaterialKind.CATALYST)
    exoi = g.input("exoi", "Exonuclease I", MaterialKind.CATALYST)
    filled = g.product(
        "filled", "nuclei with gap-filled ATAC and RNA templates",
        MaterialKind.MOLECULAR, state=both_tagged.state)
    oligo_waste = g.product("oligo_waste", "digested free single-stranded oligos",
                            MaterialKind.WASTE)
    g.react("gap_fill", "Gap-fill in nuclei with Maxima H Minus RT plus Exonuclease I",
            StepKind.EXTENSION,
            ((both_tagged, InputRole.TEMPLATE), (fill_rt, InputRole.CATALYST),
             (exoi, InputRole.CATALYST)), (filled, oligo_waste))
    g.discard(oligo_waste, "free oligos digested by Exonuclease I and removed with the supernatant")

    sorter = g.input("sort_condition",
                     "DAPI FACS into lysis wells pre-loaded with Truseq_S5_short and one N7xxx i7 primer",
                     MaterialKind.CONDITION)
    well = g.product("well", "one sorted nucleus per lysis well",
                     MaterialKind.MOLECULAR, state=filled.state)
    unsorted = g.product("unsorted", "unsorted nuclei and debris", MaterialKind.WASTE)
    g.react("sort", "Sort one DAPI-positive nucleus into each i7-barcoded lysis well",
            StepKind.PARTITION,
            ((filled, InputRole.TEMPLATE), (sorter, InputRole.CONDITION)),
            (well, unsorted), changes_molecule=False)
    g.discard(unsorted, "not deposited into a barcoded well (P24 is a no-cell control)")

    truseq_short = g.input("truseq_short", "Truseq_S5_short primer", MaterialKind.OLIGO,
                           state=short_primer("Truseq_S5_short"))
    n7 = g.input("n7_primer", "N7xxx i7 well-index primer (family)", MaterialKind.OLIGO,
                 state=n7_primer())
    q5 = g.input("q5", "NEB Q5 High-Fidelity 2x Master Mix", MaterialKind.CATALYST)
    rna_pre = g.product(
        "rna_pre", "well with i7-barcoded RNA pre-amplicons",
        MaterialKind.MOLECULAR,
        state=[*duplex_rows(ATAC_TAGGED, label="ATAC template (not yet amplified)"),
               *duplex_rows(RNA_PRE, label="i7-barcoded RNA pre-amplicon")])
    g.react("pre_pcr_rna",
            "Lyse (65 C, SDS; Tween-20 quench), gap-fill 72 C 5 min and run 4 RNA pre-PCR cycles",
            StepKind.PCR,
            ((well, InputRole.TEMPLATE), (truseq_short, InputRole.STRUCTURAL),
             (n7, InputRole.STRUCTURAL), (q5, InputRole.CATALYST)), (rna_pre,),
            geometry=PCRGeometry.exact(
                RNA_PRE, PCRPrimer("Truseq_S5_short", TRUSEQ_R1_SHORT),
                PCRPrimer("N7xxx", P7_PRIMER, anneal_nt=len(nx.S7)),
                statement="Supplementary Table 4 prints Truseq_S5_short and the N7xxx family; their 3-prime footprints face inward on the RNA product."))
    nextra_short = g.input("nextra_short", "Nextra_S5_short primer", MaterialKind.OLIGO,
                           state=short_primer("Nextra_S5_short"))
    pre = g.product(
        "pre", "well with i7-barcoded ATAC and RNA pre-amplicons",
        MaterialKind.MOLECULAR,
        state=[*duplex_rows(ATAC_PRE, label="i7-barcoded ATAC pre-amplicon"),
               *duplex_rows(RNA_PRE, label="i7-barcoded RNA pre-amplicon")])
    g.react("pre_pcr_both", "Add Nextra_S5_short and run 6 ATAC+RNA pre-PCR cycles",
            StepKind.PCR,
            ((rna_pre, InputRole.TEMPLATE), (nextra_short, InputRole.STRUCTURAL),
             (truseq_short, InputRole.STRUCTURAL), (n7, InputRole.STRUCTURAL),
             (q5, InputRole.CATALYST)), (pre,),
            geometry=PCRGeometry.exact(
                ATAC_PRE, PCRPrimer("Nextra_S5_short", NEXTRA_S5_SHORT,
                                    anneal_nt=len(NEXTRA_S5_SHORT)),
                PCRPrimer("N7xxx", P7_PRIMER, anneal_nt=len(nx.S7)),
                statement="Nextra_S5_short and N7xxx face inward on the ATAC product; the RNA pair from the previous step keeps amplifying in the same well."))

    pooled = g.product("pooled", "one plate's pooled i7-barcoded pre-amplicons",
                       MaterialKind.MOLECULAR, state=pre.state)
    column_waste = g.product("column_waste", "Zymo column flow-through",
                             MaterialKind.WASTE)
    g.react("pool", "Pool all wells of one plate and purify (Zymo DNA Clean & Concentrator)",
            StepKind.CLEANUP, ((pre, InputRole.TEMPLATE),), (pooled, column_waste),
            changes_molecule=False)
    g.discard(column_waste, "column flow-through and washes")
    exoi_post = g.input("exoi_post", "Exonuclease I (NEB 3.1 buffer)", MaterialKind.CATALYST)
    cleaned = g.product("cleaned", "pre-amplicons after residual-primer digestion",
                        MaterialKind.MOLECULAR, state=pre.state)
    primer_waste = g.product("primer_waste", "digested residual pre-PCR primers",
                             MaterialKind.WASTE)
    g.react("remove_primers",
            "Digest residual pre-PCR primers with Exonuclease I, heat-inactivate, 1.2x beads",
            StepKind.CLEANUP,
            ((pooled, InputRole.TEMPLATE), (exoi_post, InputRole.CATALYST)),
            (cleaned, primer_waste))
    g.discard(primer_waste, "digested and removed by the 1.2x bead cleanup")
    g.input("p7_primer", "illumina P7 primer", MaterialKind.OLIGO,
            state=Construct([seg("illumina P7", ILLUMINA_P7_PRIMER, "p7")],
                            name="illumina P7"))

    reviewed = [
        "nuclei", "atac_transposome", "atac_tn5", "tagment_atac",
        "capture_primer", "rt", "reverse_transcribe", "rna_transposome",
        "rna_tn5", "tagment_rna", "fill_rt", "exoi", "filled", "oligo_waste",
        "gap_fill", "sort_condition", "well", "unsorted", "sort", "truseq_short",
        "n7_primer", "q5", "rna_pre", "pre_pcr_rna", "nextra_short", "pre",
        "pre_pcr_both", "pooled", "column_waste", "pool", "exoi_post", "cleaned", "primer_waste",
        "remove_primers", "p7_primer",
    ]
    # Steps 61-63: one 8-cycle tube with ATAC_plate_S5xx, RNA_plate_S5xx and Illumina
    # P7.  Drawn as two primer-pair readings of the same tube so each product keeps
    # exact geometry; there is no physical split.
    atac_geometry = PCRGeometry.exact(
        ATAC, PCRPrimer("ATAC_plate_SS5xx", ATAC_P5_PRIMER, anneal_nt=len(nx.S5)),
        PCRPrimer("illumina P7", ILLUMINA_P7_PRIMER),
        statement="Supplementary Table 4 defines ATAC_plate_SS5xx and the illumina P7 primer; they face inward on the i7-barcoded ATAC pre-amplicon.")
    rna_geometry = PCRGeometry.exact(
        RNA, PCRPrimer("RNA_plate_S5xx", RNA_P5_PRIMER, anneal_nt=18),
        PCRPrimer("illumina P7", ILLUMINA_P7_PRIMER),
        statement="Supplementary Table 4 defines RNA_plate_S5xx (3-prime 18 nt within the capture arm) and the illumina P7 primer; they face inward on the RNA pre-amplicon.")
    for prefix, label, final, primers, geometry in (
            ("atac", "plate ISSAAC-seq ATAC", ATAC, ATAC_PRIMERS, atac_geometry),
            ("rna", "plate ISSAAC-seq RNA", RNA, RNA_PRIMERS, rna_geometry)):
        forward = g.input(prefix + "_forward",
                          ("ATAC_plate_SS5xx" if prefix == "atac" else "RNA_plate_S5xx")
                          + " i5 plate/modality primer", MaterialKind.OLIGO,
                          state=plate_i5_primer(prefix))
        library = g.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
        g.react(prefix + "_pcr",
                "8-cycle final PCR (one tube with both i5 primers): " + label + " product",
                StepKind.PCR,
                ((cleaned, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 ("p7_primer", InputRole.STRUCTURAL), ("q5", InputRole.CATALYST)),
                (library,), geometry=geometry)
        g.sequence(library, SequencingPlan(
            library=final, primers=tuple(primers),
            acquired_roles=tuple(p.role for p in primers),
            read_lengths=PLATE_READ_LENGTHS,
            read_uses=READ_USES[("Plate ATAC library" if prefix == "atac"
                                 else "Plate RNA library")]))
        reviewed.extend((forward.id, prefix + "_pcr", library.id))
    g.declare_reviewed_scope(reviewed, SOURCE)
    return g.require_publishable()
