"""INVADEseq: targeted bacterial 16S capture beside 10x 5' host RNA capture."""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp, oligo

TITLE = "INVADEseq — invasion–adhesion-directed expression sequencing"
NOTES = "01_invadeseq.html"
SOURCE = ('Defining protocol: <a href="https://doi.org/10.1038/s41596-023-00888-7">'
          'Galeano Niño et al., <i>Nature Protocols</i> (2023)</a>.')
SUMMARY = ("A bacterial 16S reverse-transcription primer is spiked into a 10x Chromium "
           "5′ single-cell reaction beside the ordinary poly(dT) primer. Host and bacterial "
           "cDNAs therefore inherit the same cell barcode; an aliquot is then amplified "
           "twice with nested 16S primers to make a separate bacterial library.")
CAVEAT = ("The three INVADEseq primers are published exactly. The underlying 10x gel-bead "
          "reagent is commercial: its published barcode, UMI and adapter roles are shown, "
          "but undisclosed bead chemistry is not reconstructed.")

RT_16S = "GGGTTGCGCTCGTTG"
ENRICH_FORWARD = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC"
ENRICH_REVERSE = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTTCACGRCACGAGCTGACGAC"
ENRICH_FORWARD_BINDING = ENRICH_FORWARD[len(il.P5):]
ENRICH_REVERSE_BINDING = ENRICH_REVERSE[len(il.TRUSEQ_READ2):]

CELL_BARCODE = feature("cell_barcode", "cell_barcode", "whitelist",
                       whitelist="10x Chromium 5-prime v1/v2 whitelist")
UMI = feature("umi", "umi", "random")
I5 = feature("sample_index_i5", "sample_index", "unknown")
I7 = feature("sample_index_i7", "sample_index", "unknown")
TSO = "TTTCTTATATGGG"


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name, top, tag, **kw)


def oligos():
    return (
        oligo("16S 1100R RT primer", [seg("16S 1100R", RT_16S)]),
        oligo("16S Enrich Forward", [seg("P5", il.P5, "p5"),
                                      seg("10x-end binding region", ENRICH_FORWARD_BINDING, "r1")]),
        oligo("1061R Enrich Nested Reverse", [seg("Read 2 arm", il.TRUSEQ_READ2, "r2"),
                                               seg("nested 1061R", ENRICH_REVERSE_BINDING,
                                                   placeholder=True, length_bp=20,
                                                   note="published IUPAC R is A or G")]),
    )


def input_scene() -> Scene:
    sc = Scene()
    sc.strand("host", [seg("host mRNA", "X" * 28, placeholder=True),
                       seg("poly(A)", "A" * 18)], label="host RNA")
    sc.strand("bacterial", [seg("bacterial 16S rRNA", "X" * 34, placeholder=True),
                            seg("1100R target", revcomp(RT_16S))], label="bacterial RNA")
    return sc


def dual_rt_scene() -> Scene:
    sc = input_scene()
    sc.anneal("poly(dT) primer", [seg("poly(dT)", "T" * 18)], to="host",
              pair=("poly(dT)", "poly(A)"), label="standard poly(dT) RT primer")
    sc.arrow("poly(dT) primer", "reverse transcriptase")
    sc.anneal("16S primer", [seg("16S 1100R", RT_16S)], to="bacterial",
              pair=("16S 1100R", "1100R target"), label="16S 1100R")
    sc.arrow("16S primer", "reverse transcriptase")
    return sc


def barcoded_cdna_scene() -> Scene:
    common = [seg("cell barcode", "B" * 16, "cbc", placeholder=True,
                  feature=CELL_BARCODE),
              seg("UMI", "U" * 10, "umi", placeholder=True, feature=UMI),
              seg("template-switch sequence", TSO, "tso")]
    sc = Scene()
    sc.strand("host", [*common, seg("host cDNA", "X" * 36, placeholder=True)],
              label="host cDNA")
    sc.strand("bacterial", [*common, seg("16S cDNA", "X" * 36, placeholder=True)],
              label="bacterial cDNA")
    sc.labels("host"); sc.labels("bacterial")
    return sc


def split_scene() -> Scene:
    sc = Scene()
    sc.strand("gex", [seg("GEX aliquot", "X" * 30, placeholder=True)],
              label="standard 5′ GEX branch")
    sc.strand("16s", [seg("16S aliquot", "X" * 30, placeholder=True)],
              label="bacterial enrichment branch")
    sc.note("16s", "two nested PCRs; retain 955–1,215 bp")
    return sc


ENRICHED_AMPLICON = Construct([
    seg("P5 from 16S Enrich Forward", il.P5, "p5"),
    seg("10x-end binding region", ENRICH_FORWARD_BINDING, "r1"),
    seg("cell barcode", "B" * 16, "cbc", placeholder=True, feature=CELL_BARCODE),
    seg("UMI", "U" * 10, "umi", placeholder=True, feature=UMI),
    seg("template-switch sequence", TSO, "tso"),
    seg("16S cDNA to nested 1061R", "X" * 52, placeholder=True),
    seg("reverse-complemented nested-primer end", "X" * 20, placeholder=True,
        length_bp=20),
    seg("PCR-completed Read 2 leading base", "A", "r2"),
    seg("Read 2 arm", il.INDEX1_PRIMER, "r2"),
], name="twice-enriched INVADEseq 16S amplicon")


def enriched_scene() -> Scene:
    sc = Scene.duplex(list(ENRICHED_AMPLICON), label="16S enrichment amplicon")
    sc.strand("gex", [seg("reserved standard 5-prime GEX aliquot", "X" * 32,
                          placeholder=True)], label="host GEX branch")
    sc.note("top", "two PCR rounds with the published enrichment primers")
    sc.note("bottom", "BluePippin retains complete 955–1,215-bp products")
    return sc


def indexed_scene() -> Scene:
    sc = Scene()
    sc.strand("gex", [seg("dual-index host GEX library", "X" * 38,
                          placeholder=True)], label="host GEX library")
    sc.strand("16s", [seg("dual-index 16S enrichment library", "X" * 38,
                          placeholder=True)], label="bacterial 16S library")
    sc.note("gex", "high-throughput Illumina sequencing")
    sc.note("16s", "2 × 300 MiSeq v3")
    return sc


def common_left() -> list[Segment]:
    return [seg("P5", il.P5, "p5"),
            seg("i5", "J" * 10, "cbc", placeholder=True, feature=I5),
            seg("Read 1 arm", il.TRUSEQ_READ1, "r1"),
            seg("cell barcode", "B" * 16, "cbc", placeholder=True, feature=CELL_BARCODE),
            seg("UMI", "U" * 10, "umi", placeholder=True, feature=UMI),
            seg("template-switch sequence", TSO, "tso")]


def indexed_right() -> list[Segment]:
    return [seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
            seg("i7 reverse complement", "I" * 10, "cbc", placeholder=True, feature=I7),
            seg("P7 reverse complement", il.P7_RC, "p7")]


GEX_LIBRARY = Construct([
    *common_left(), seg("host cDNA fragment", "X" * 42, placeholder=True),
    seg("dA ligation junction", "A"), *indexed_right()], name="INVADEseq host GEX library")

BACTERIAL_LIBRARY = Construct([
    *common_left(),
    seg("16S cDNA through nested 1061R site", "X" * 64, placeholder=True,
        note="length varies with the captured 16S molecule; BluePippin selects 955–1,215-bp products"),
    # The full published nested reverse primer supplies TruSeq Read 2.  On this displayed
    # strand its reverse complement is one leading A plus INDEX1_PRIMER; unlike the GEX
    # branch, that A is PCR-installed and is not a ligation junction.
    seg("PCR-completed Read 2 leading base", "A", "r2"),
    *indexed_right()], name="INVADEseq bacterial 16S enrichment library")

PRIMERS = tuple(sp.TRUSEQ[k] for k in ("R1", "I1", "I2", "R2"))
for _name, _library in (("GEX", GEX_LIBRARY), ("bacterial 16S", BACTERIAL_LIBRARY)):
    _problems = sp.verify(_library, PRIMERS, required_roles=tuple(p.role for p in PRIMERS))
    if _problems:
        raise ValueError(f"INVADEseq {_name} library: " + "; ".join(_problems))

FINAL_LIBRARIES = (
    ("Final host gene-expression library", GEX_LIBRARY, PRIMERS,
     "The ordinary 10x 5′ GEX branch reports host transcripts with the same cell barcode space used by the bacterial branch.",
     "Standard TruSeq sequencing primers bind the completed dual-index 10x library."),
    ("Final bacterial 16S enrichment library", BACTERIAL_LIBRARY, PRIMERS,
     "Two targeted PCRs retain the 10x cell barcode and UMI while enriching cDNA between the bead end and nested 1061R site.",
     "Read 1 begins with the 16-base cell barcode and 10-base UMI before continuing through the template-switch sequence into 16S cDNA; Read 2 enters 16S cDNA from the nested-primer end."),
)

READ_LENGTHS = {
    "Final host gene-expression library": {
        "Read 1": 26, "Index 1 (i7)": 10, "Index 2 (i5)": 10, "Read 2": 90},
    "Final bacterial 16S enrichment library": {
        "Read 1": 300, "Index 1 (i7)": 10, "Index 2 (i5)": 10, "Read 2": 300},
}


def sections():
    return [
        ("Co-encapsulate one host cell and its associated bacteria",
         input_scene().rows(),
         "The host cell defines the single-cell partition; associated bacterial RNA enters the same GEM."),
        ("Prime host poly(A) RNA and bacterial 16S rRNA in the same GEM",
         dual_rt_scene().rows(),
         "The published 16S 1100R spike-in is added beside the standard 10x poly(dT) RT primer."),
        ("Template-switch both products onto the barcoded gel-bead oligo",
         barcoded_cdna_scene().rows(),
         "The exact commercial bead chemistry is undisclosed; the source nevertheless establishes that both products acquire the same 10x cell barcode and independent UMIs."),
        ("Amplify cDNA, then split host GEX and bacterial-enrichment branches",
         split_scene().rows(),
         "The GEX aliquot follows standard 10x library preparation; the second aliquot enters targeted enrichment."),
        ("Enrich the barcode-bearing 16S molecule twice and size-select it",
         enriched_scene().rows(),
         "The bacterial aliquot undergoes 35-cycle enrichment, 955–1,215-bp selection, then a second 20-cycle enrichment with the same published primer pair."),
        ("Add dual sample indexes to both completed libraries",
         indexed_scene().rows(),
         "Dual Index TT Set A PCR completes both Illumina libraries; the bacterial branch is size-selected again before 2 × 300 MiSeq sequencing."),
    ]
