"""DR-seq co-amplification and source-backed branch intermediates."""
from __future__ import annotations

from adapter_ligation import AdapterLigationGeometry
from batch_ngs import seg
from chemdraw import (Construct, MolecularState, Row, Scene, Workflow,
                      complement_segments, duplex_rows, feature)
from multibranch import illumina_branch
from multimodal_spatial import finish_pcr_branch
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind

TITLE = "DR-seq — genome and transcriptome without physical separation"
NOTES = "01_dr-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/nbt.3129">'
          'Dey et al., <i>Nature Biotechnology</i> (2015)</a>.')
SUMMARY = ("A cell-barcoded oligo(dT) first makes cDNA in the same lysate as genomic DNA. "
           "Ad-2 random primers then quasilinearly copy both analytes before the pool is "
           "split into genome and transcriptome finishing reactions.")
CAVEAT = ("The paper prints Ad-1x, Ad-2 and the branch-selection primers, but refers the "
          "RNA branch to the original CEL-Seq small-RNA finishing workflow without printing "
          "those outer library oligos. The final shells therefore show inferred standard "
          "Illumina adapters; the co-amplification intermediates are source-exact.")

AD1_SPACER = "CGATTGAGGCCGG"
T7_PROMOTER = "TAATACGACTCACTATAGGG"
RA5 = "GTTCAGAGTTCTACAGTCCGACGATC"
CELL_BARCODE_NT = 8
POLYT_NT = 24
AD2_HANDLE = "GTGAGTGATGGTTGAGGTAGTGTGGAG"
AD2_RANDOM_NT = 8
P1 = "CGATTGAGGCCGGTAATAC"
P2 = AD2_HANDLE
P3 = "GTGAGCTGGAGTTGAGGTAGTGTGGAG"


def ad1_segments():
    return [seg("Ad-1x spacer", AD1_SPACER), seg("T7 promoter", T7_PROMOTER, "t7"),
            seg("RA5", RA5, "r1"),
            seg("cell barcode", "B" * CELL_BARCODE_NT, "cbc", placeholder=True,
                feature=feature("cell_barcode", "cell_barcode", "fixed")),
            seg("oligo(dT)24", "T" * POLYT_NT),
            seg("V anchor", "V", placeholder=True)]


def ad2_segments():
    return [seg("Ad-2 common handle", AD2_HANDLE, "tso"),
            seg("Ad-2 random N8", "N" * AD2_RANDOM_NT, placeholder=True)]


def rt_scene() -> Scene:
    mrna = [seg("transcript", "X" * 34, placeholder=True),
            seg("poly(A)", "A" * POLYT_NT)]
    sc = Scene(); sc.strand("mRNA", mrna, label="poly(A) RNA")
    sc.anneal("Ad-1x", ad1_segments(), to="mRNA",
              pair=("oligo(dT)24", "poly(A)"), label="cell-barcoded Ad-1x",
              unpaired=("V anchor",))
    sc.arrow("Ad-1x", "reverse transcription")
    sc.mark("Ad-1x", "cell barcode", "8-nt cell identity")
    sc.labels("Ad-1x")
    return sc


def ad2_priming_scene(template_name: str) -> Scene:
    template = [seg(template_name, "X" * 40, placeholder=True),
                seg("random priming site", "NNNNNNNN", placeholder=True)]
    sc = Scene(); sc.strand("template", template, label=template_name)
    sc.anneal("Ad-2", ad2_segments(), to="template",
              pair=("Ad-2 random N8", "random priming site"), label="Ad-2 random primer")
    sc.arrow("Ad-2", "strand-displacing extension")
    sc.labels("Ad-2")
    return sc


def gdna_amplicon() -> Construct:
    left = [*ad2_segments(), seg("genomic copy", "X" * 40, placeholder=True)]
    return Construct([*left, *complement_segments(ad2_segments(), suffix=" right complement")],
                     name="Ad-2/Ad-2 genomic full amplicon")


def cdna_amplicon() -> Construct:
    # The productive transcript-derived species is Ad-2 at its random-primed end and
    # retains Ad-1x at the poly(A)-primed end.
    return Construct([*ad2_segments(), seg("cDNA", "X" * 38, placeholder=True),
                      seg("copied poly(A)", "A" * POLYT_NT),
                      seg("cell barcode complement", "b" * CELL_BARCODE_NT, "cbc",
                          placeholder=True,
                          feature=feature("cell_barcode", "cell_barcode", "fixed")),
                      *complement_segments(ad1_segments()[:3], suffix=" complement")],
                     name="Ad-2/Ad-1x transcript amplicon")


GDNA_AMPLICON = gdna_amplicon()
CDNA_AMPLICON = cdna_amplicon()


def arna_product() -> Construct:
    return Construct([
        seg("T7 initial G-rich end", "GGG", "t7"), seg("RA5", RA5, "r1"),
        seg("cell barcode", "B" * CELL_BARCODE_NT, "cbc", placeholder=True,
            feature=feature("cell_barcode", "cell_barcode", "fixed")),
        seg("poly(U) and antisense transcript", "UUU" + "X" * 38, placeholder=True),
        *complement_segments(ad2_segments(), suffix=" complement"),
    ], name="cell-barcoded antisense RNA")


ARNA = arna_product()


def arna_rows():
    sc = Scene(); sc.strand("aRNA", list(ARNA), label="cell-barcoded aRNA")
    sc.labels("aRNA")
    return sc.rows()

# The exact post-sonication DNA insert and exact post-IVT RNA insert are retained, while
# outer Illumina completion is marked inferred because those oligos are not printed here.
DNA, DP = illumina_branch("DR-seq genome", "sonicated genomic insert",
                          architecture="truseq", inferred_adapters=True)
RNA, RP = illumina_branch("DR-seq transcriptome", "3-prime RNA-derived insert",
                          cell=CELL_BARCODE_NT, inferred_adapters=True,
                          encoding="fixed")
FINAL_LIBRARIES = (
    ("Genome library", DNA, DP,
     "NEBNext Ultra library from sonicated, Ad-2-depleted genomic amplicons; outer shell is inferred because its lot-specific oligos are not printed.",
     "Canonical TruSeq-family primer landing on the inferred NEBNext shell."),
    ("Transcriptome library", RNA, RP,
     "Directional CEL-Seq-derived library retaining the 8-nt Ad-1x cell barcode; outer small-RNA completion is shown as inferred.",
     "Primer landing is shown on the inferred Illumina shell; see the evidence caveat."),
)


def oligos():
    return (
        f"<b>Ad-1x</b> — 5′-{AD1_SPACER}{T7_PROMOTER}{RA5}[8-nt cell barcode]{'T' * POLYT_NT}V-3′",
        f"<b>Ad-2</b> — 5′-{AD2_HANDLE}[8 random nt]-3′",
        f"<b>P1</b> — 5′-{P1}-3′",
        f"<b>P2</b> — 5′-{P2}-3′",
        f"<b>5′-biotin P3</b> — 5′-biotin-{P3}-3′",
    )


def workflow() -> Workflow:
    start = Scene.duplex([seg("genomic DNA", "X" * 42, placeholder=True)],
                         label="single-cell genomic DNA")
    wf = Workflow(MolecularState("One lysed cell: genomic DNA plus poly(A) RNA",
                                 tuple(start.rows())))
    wf.react("Reverse-transcribe RNA with cell-barcoded Ad-1x",
             rt_scene().rows(), name="Ad-1x-tagged cDNA beside genomic DNA",
             note="The 5′ tail installs the T7 promoter, RA5 and cell barcode; genomic DNA remains in the tube.")
    wf.react("Random-prime genomic DNA and cDNA with Ad-2",
             [*ad2_priming_scene("genomic DNA").rows(),
              *ad2_priming_scene("Ad-1x-tagged cDNA").rows()],
             name="First Ad-2 semi-amplicons")
    wf.react("Run seven MALBAC-like quasilinear cycles",
             [*duplex_rows(GDNA_AMPLICON, label="genome-derived Ad-2/Ad-2 species"),
              *duplex_rows(CDNA_AMPLICON, label="transcript-derived Ad-2/Ad-1x species"),
              Row(chunks=[("full Ad-2/Ad-2 amplicons self-anneal through complementary ends and leave the quasilinear pool", None, False)])],
             name="Co-amplified genome and transcript products")
    wf.react("Split the common amplification pool into DNA and RNA aliquots",
             [Row(chunks=[("DNA aliquot → P2 PCR → 5′-biotin P3 PCR → sonicate → streptavidin depletion", None, False)]),
              Row(chunks=[("RNA aliquot → P1 second strand → T7 IVT", None, False)])],
             name="Modality-specific finishing reactions")
    wf.react("Select the RNA branch with P1 and amplify it by T7 transcription",
             [*duplex_rows(CDNA_AMPLICON, label="P1-completed T7 template"),
              Row(chunks=[("T7 RNA polymerase →", "t7", False)]),
              *arna_rows()],
             name="Antisense RNA copies from transcript-derived templates",
             note="Genomic Ad-2/Ad-2 products lack the Ad-1x T7 end and do not make this aRNA.")
    wf.react("Finish the DNA and RNA aliquots as separate Illumina libraries",
             [*duplex_rows(DNA, label="genome library"),
              *duplex_rows(RNA, label="transcriptome library")],
             name="DR-seq genome and transcriptome libraries",
             note="The source names NEBNext Ultra for DNA and the CEL-Seq workflow for RNA; the unprinted outer oligos are marked inferred.")
    return wf


def protocol_graph():
    """Typed provenance for DR-seq's common amplification and two finishes."""
    g = ProtocolGraph(TITLE, authority=SOURCE)
    cell = g.input(
        "cell_lysate", "one lysed cell containing genomic DNA and poly(A) RNA",
        MaterialKind.BIOLOGICAL,
        state=[*Scene.duplex([seg("genomic DNA", "X" * 42, placeholder=True)],
                             label="single-cell genomic DNA").rows(),
               *rt_scene().rows()], primary=True)

    ad1 = g.input("ad1", "cell-barcoded Ad-1x oligo(dT)24V primer",
                  MaterialKind.OLIGO, state=rt_scene())
    rt = g.input("rt", "reverse transcriptase", MaterialKind.CATALYST)
    copied = g.product(
        "copied", "Ad-1x-tagged cDNA beside genomic DNA",
        MaterialKind.MOLECULAR,
        state=[*rt_scene().rows(),
               *Scene.duplex([seg("genomic DNA", "X" * 42, placeholder=True)],
                             label="unchanged genomic DNA").rows()])
    g.react(
        "reverse_transcribe", "Reverse-transcribe RNA with cell-barcoded Ad-1x",
        StepKind.REVERSE_TRANSCRIPTION,
        ((cell, InputRole.TEMPLATE), (ad1, InputRole.STRUCTURAL),
         (rt, InputRole.CATALYST)), (copied,))

    ad2 = g.input("ad2", "Ad-2 common-handle random-N8 primer",
                  MaterialKind.OLIGO)
    strand_polymerase = g.input(
        "strand_polymerase", "strand-displacing DNA polymerase",
        MaterialKind.CATALYST)
    primed = g.product(
        "ad2_primed", "Ad-2-primed genomic DNA and cDNA copies",
        MaterialKind.MOLECULAR,
        state=[*ad2_priming_scene("genomic DNA").rows(),
               *ad2_priming_scene("Ad-1x-tagged cDNA").rows()])
    g.react(
        "random_prime", "Random-prime genomic DNA and cDNA with Ad-2",
        StepKind.EXTENSION,
        ((copied, InputRole.TEMPLATE), (ad2, InputRole.STRUCTURAL),
         (strand_polymerase, InputRole.CATALYST)), (primed,))
    coamplified = g.product(
        "coamplified", "quasilinearly co-amplified genome and transcript products",
        MaterialKind.MOLECULAR,
        state=[*duplex_rows(GDNA_AMPLICON, label="genome-derived Ad-2/Ad-2 species"),
               *duplex_rows(CDNA_AMPLICON,
                            label="transcript-derived Ad-2/Ad-1x species")])
    g.react(
        "quasilinear_cycles", "Run seven MALBAC-like quasilinear extension cycles",
        StepKind.EXTENSION,
        ((primed, InputRole.TEMPLATE), (ad2, InputRole.STRUCTURAL),
         (strand_polymerase, InputRole.CATALYST)), (coamplified,))

    dna = g.product("dna_branch", "genome finishing aliquot",
                    MaterialKind.MOLECULAR, state=GDNA_AMPLICON)
    rna = g.product("rna_branch", "transcriptome finishing aliquot",
                    MaterialKind.MOLECULAR, state=CDNA_AMPLICON)
    g.react("split", "Split the common amplification pool into DNA and RNA aliquots",
            StepKind.SPLIT, ((coamplified, InputRole.TEMPLATE),), (dna, rna),
            changes_molecule=False)

    p2_forward = g.input("p2_forward", "P2 genome-amplification primer",
                         MaterialKind.OLIGO)
    p2_reverse = g.input("p2_reverse", "opposing P2 genome-amplification primer",
                         MaterialKind.OLIGO)
    p2_polymerase = g.input("p2_polymerase", "P2 PCR polymerase",
                            MaterialKind.CATALYST)
    p2_product = g.product(
        "p2_product", "P2-amplified genomic products", MaterialKind.MOLECULAR,
        state=GDNA_AMPLICON)
    g.react(
        "p2_pcr", "Amplify the genome aliquot with P2", StepKind.PCR,
        ((dna, InputRole.TEMPLATE), (p2_forward, InputRole.STRUCTURAL),
         (p2_reverse, InputRole.STRUCTURAL),
         (p2_polymerase, InputRole.CATALYST)), (p2_product,),
        geometry=PCRGeometry.unavailable(
            GDNA_AMPLICON,
            reason="The source gives one P2 sequence for the repeated Ad-2 handles; this representative construct does not assert a unique oriented primer pair."))

    p3 = g.input("p3_biotin", "5-prime-biotin P3 primer", MaterialKind.OLIGO)
    p3_partner = g.input("p3_partner", "opposing genome PCR primer",
                         MaterialKind.OLIGO, inferred=True,
                         source=SOURCE + " (opposing primer identity not resolved in this page model)")
    p3_product = g.product(
        "p3_product", "biotin-P3-amplified genomic products",
        MaterialKind.MOLECULAR, state=GDNA_AMPLICON, inferred=True,
        source=SOURCE + " (opposing primer identity not resolved in this page model)")
    g.react(
        "p3_pcr", "Amplify the genome branch with 5-prime-biotin P3",
        StepKind.PCR,
        ((p2_product, InputRole.TEMPLATE), (p3, InputRole.STRUCTURAL),
         (p3_partner, InputRole.STRUCTURAL),
         (p2_polymerase, InputRole.CATALYST)), (p3_product,),
        geometry=PCRGeometry.unavailable(
            GDNA_AMPLICON,
            reason="The source prints biotin-P3, but this page has not resolved a unique opposing primer landing site."))
    shear = g.input("shearing", "sonication condition", MaterialKind.CONDITION)
    fragmented = g.product("fragmented_dna", "sonicated genome amplicons",
                           MaterialKind.MOLECULAR, state=p3_product.state)
    g.react("sonicate", "Sonicate the biotin-labelled genome amplicons",
            StepKind.FRAGMENTATION,
            ((p3_product, InputRole.TEMPLATE), (shear, InputRole.CONDITION)),
            (fragmented,))
    streptavidin = g.input("streptavidin", "streptavidin depletion matrix",
                           MaterialKind.SELECTION)
    depleted = g.product("depleted_dna", "Ad-2-depleted sonicated DNA",
                         MaterialKind.MOLECULAR, state=fragmented.state)
    bound = g.product("bound_ad2", "streptavidin-bound Ad-2 end fragments",
                      MaterialKind.WASTE)
    g.react(
        "deplete_ad2", "Deplete biotin-labelled Ad-2 ends with streptavidin",
        StepKind.SELECTION,
        ((fragmented, InputRole.TEMPLATE),
         (streptavidin, InputRole.SELECTION)), (depleted, bound),
        changes_molecule=False)
    g.discard(bound, "contains the biotin-labelled Ad-2 ends removed before library prep")
    neb_adapter = g.input(
        "neb_adapter", "NEBNext Ultra DNA-library adapter",
        MaterialKind.ADAPTER, inferred=True,
        source=SOURCE + " (kit named; lot-specific adapter oligos not printed)")
    dna_adapted = g.product(
        "dna_adapted", "NEBNext-adapted genome fragments",
        MaterialKind.MOLECULAR, state=DNA, inferred=True,
        source=SOURCE + " (outer NEBNext shell inferred from the named kit)")
    g.react(
        "adapt_dna", "Ligate the named NEBNext Ultra adapter",
        StepKind.LIGATION,
        ((depleted, InputRole.TEMPLATE), (neb_adapter, InputRole.STRUCTURAL)),
        (dna_adapted,), geometry=AdapterLigationGeometry.unavailable(
            "The paper names NEBNext Ultra but does not print its lot-specific adapter-end geometry."))

    p1 = g.input("p1", "P1 transcript-branch primer", MaterialKind.OLIGO)
    second_polymerase = g.input("second_polymerase", "second-strand polymerase",
                                MaterialKind.CATALYST)
    t7_template = g.product(
        "t7_template", "P1-completed cell-barcoded T7 template",
        MaterialKind.MOLECULAR, state=CDNA_AMPLICON)
    g.react(
        "second_strand", "Select transcript-derived products with P1 and make second strand",
        StepKind.EXTENSION,
        ((rna, InputRole.TEMPLATE), (p1, InputRole.STRUCTURAL),
         (second_polymerase, InputRole.CATALYST)), (t7_template,))
    t7_polymerase = g.input("t7_polymerase", "T7 RNA polymerase",
                            MaterialKind.CATALYST)
    arna = g.product("arna", "cell-barcoded antisense RNA copies",
                     MaterialKind.MOLECULAR, state=ARNA)
    g.react("ivt", "Amplify the transcript branch by T7 transcription",
            StepKind.TRANSCRIPTION,
            ((t7_template, InputRole.TEMPLATE),
             (t7_polymerase, InputRole.CATALYST)), (arna,))
    cel_adapter = g.input(
        "cel_adapter", "CEL-Seq small-RNA finishing adapter",
        MaterialKind.ADAPTER, inferred=True,
        source=SOURCE + " (workflow referenced; outer oligos not printed)")
    rna_adapted = g.product(
        "rna_adapted", "CEL-Seq-finished antisense-RNA library template",
        MaterialKind.MOLECULAR, state=RNA, inferred=True,
        source=SOURCE + " (outer small-RNA shell inferred from the referenced workflow)")
    g.react(
        "adapt_rna", "Attach the referenced CEL-Seq small-RNA finishing adapter",
        StepKind.LIGATION,
        ((arna, InputRole.TEMPLATE), (cel_adapter, InputRole.STRUCTURAL)),
        (rna_adapted,), geometry=AdapterLigationGeometry.unavailable(
            "The defining paper refers to the CEL-Seq finishing workflow without printing its outer adapter oligos."))

    reviewed = [
        "cell_lysate", "ad1", "rt", "reverse_transcribe", "ad2",
        "strand_polymerase", "random_prime", "quasilinear_cycles", "split",
        "p2_forward", "p2_reverse", "p2_polymerase", "p2_pcr", "p3_biotin",
        "p3_partner", "p3_pcr", "shearing", "sonicate", "streptavidin",
        "deplete_ad2", "neb_adapter", "adapt_dna", "p1", "second_polymerase",
        "second_strand", "t7_polymerase", "ivt", "cel_adapter", "adapt_rna",
    ]
    finish_pcr_branch(
        g, "dna_final", dna_adapted, DNA, DP, reviewed,
        label="DR-seq genome", inferred_primers=True,
        limitation="The defining paper names NEBNext Ultra but does not print a universal lot-specific PCR-primer pair.")
    finish_pcr_branch(
        g, "rna_final", rna_adapted, RNA, RP, reviewed,
        label="DR-seq transcriptome", inferred_primers=True,
        limitation="The defining paper refers to CEL-Seq finishing without printing the outer PCR-primer pair.")
    g.declare_reviewed_scope(reviewed, SOURCE)
    return g.require_publishable()
