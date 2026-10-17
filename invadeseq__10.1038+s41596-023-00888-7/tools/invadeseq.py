"""INVADEseq: targeted bacterial 16S capture beside 10x 5' host RNA capture."""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp, oligo
from multimodal_spatial import finish_pcr_branch
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "INVADEseq — invasion–adhesion-directed expression sequencing"
NOTES = "01_invadeseq.html"
SOURCE = ('Defining protocol: <a href="https://doi.org/10.1038/s41596-023-00888-7">'
          'Galeano Niño et al., <i>Nature Protocols</i> (2023)</a>.')
SUMMARY = ("A bacterial 16S reverse-transcription primer is spiked into a 10x Chromium "
           "5′ single-cell reaction beside the ordinary poly(dT) primer. Host and bacterial "
           "cDNAs therefore inherit the same cell barcode; an aliquot is then amplified "
           "twice with nested 16S primers to make a separate bacterial library.")
CAVEAT = ("The three INVADEseq primers are published exactly. The cited 10x CG000086 "
          "Rev J guide remains unverified; a later official Rev K guide prints the gel-bead "
          "primer architecture. Exact reagent-revision equivalence and selected sample "
          "index member remain unresolved. The protocol Materials list (lines 154–156) "
          "identifies v2 reagent kit PN-1000263 and library kit PN-1000002, but Step 85 "
          "(line 758) calls PN-1000002 a v2 kit. This page uses the Materials identities "
          "for those two roles; which kit combination the Step 85 wording intended remains unclear. "
          "Host GEX read lengths (26 × 10 × 10 × 90) come from the 2022 <i>Nature</i> application "
          "supplement (MOESM1 Methods, NextSeq 2000 P3-100); the 2023 protocol instead names HiSeq 4000 "
          "paired-end 75-bp reads in its Equipment setup (line 306) and leaves NextSeq 2000/NovaSeq "
          "read lengths to the user in Step 84 (line 732), so GEX cycles are application-specific, not "
          "a protocol requirement. The bacterial library is indexed with 15 cycles and then re-selected "
          "with the same 955–1,215-bp BluePippin program (Steps 115 and 119); the post-index size "
          "change is not computed here.")

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
    seg("reverse-complemented nested-primer end", revcomp(ENRICH_REVERSE_BINDING),
        placeholder=True, length_bp=20,
        note="published IUPAC R is A or G"),
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
    seg("16S cDNA before nested 1061R site", "X" * 44, placeholder=True,
        note="length varies with the captured 16S molecule; BluePippin selects 955–1,215-bp products"),
    seg("reverse-complemented nested 1061R site", revcomp(ENRICH_REVERSE_BINDING),
        placeholder=True, length_bp=20,
        note="published IUPAC R is A or G"),
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
         "Dual Index TT Set A PCR completes both Illumina libraries (14–16 cycles for GEX, 15 for the bacterial branch); the bacterial branch is size-selected again with the same 955–1,215-bp BluePippin program before 2 × 300 MiSeq v3 sequencing."),
    ]


def _enrichment_geometry() -> PCRGeometry:
    return PCRGeometry.exact(
        ENRICHED_AMPLICON,
        PCRPrimer("16S Enrich Forward", ENRICH_FORWARD),
        PCRPrimer("1061R Enrich Nested Reverse", ENRICH_REVERSE),
        statement="The protocol prints both enrichment primers; the published 3-prime 1061R footprint and the 10x-end footprint face inward on the enriched amplicon.")


def protocol_graph():
    g = ProtocolGraph(TITLE, authority=SOURCE)
    associated = g.input(
        "associated_cells", "host cells with associated bacterial RNA",
        MaterialKind.BIOLOGICAL, state=input_scene(), primary=True)
    gem_condition = g.input("gem_condition", "10x GEM partitioning condition",
                            MaterialKind.CONDITION)
    gem = g.product("gem", "one host-cell partition containing associated bacterial RNA",
                    MaterialKind.MOLECULAR, state=input_scene())
    excluded = g.product("excluded", "empty or colliding GEM material",
                         MaterialKind.WASTE)
    g.react("partition", "Co-encapsulate a host cell and its associated bacteria",
            StepKind.PARTITION,
            ((associated, InputRole.TEMPLATE),
             (gem_condition, InputRole.CONDITION)), (gem, excluded),
            changes_molecule=False)
    g.discard(excluded, "does not preserve an interpretable single-host-cell partition")

    poly_dt = g.input("poly_dt", "standard 10x poly(dT) RT primer",
                      MaterialKind.OLIGO, inferred=True,
                      source=SOURCE + " (commercial bead chemistry undisclosed)")
    primer_16s = g.input("primer_16s", "published 16S 1100R RT primer",
                         MaterialKind.OLIGO)
    rt = g.input("rt", "reverse transcriptase", MaterialKind.CATALYST)
    copied = g.product(
        "copied", "host and bacterial first-strand cDNAs in one GEM",
        MaterialKind.MOLECULAR, state=dual_rt_scene())
    g.react("reverse_transcribe",
            "Prime host poly(A) RNA and bacterial 16S rRNA in the same GEM",
            StepKind.REVERSE_TRANSCRIPTION,
            ((gem, InputRole.TEMPLATE), (poly_dt, InputRole.STRUCTURAL),
             (primer_16s, InputRole.STRUCTURAL), (rt, InputRole.CATALYST)),
            (copied,))

    bead_oligo = g.input(
        "bead_oligo", "10x cell-barcode/UMI/template-switch gel-bead oligo",
        MaterialKind.OLIGO, inferred=True,
        source=SOURCE + " (commercial bead chemistry undisclosed)")
    switched = g.product(
        "switched", "host and bacterial cDNAs carrying the same cell barcode",
        MaterialKind.MOLECULAR, state=barcoded_cdna_scene(), inferred=True,
        source=SOURCE + " (barcode roles established; bead chemistry undisclosed)")
    g.react("template_switch",
            "Template-switch both products onto the barcoded gel-bead oligo",
            StepKind.TEMPLATE_SWITCHING,
            ((copied, InputRole.TEMPLATE), (bead_oligo, InputRole.STRUCTURAL),
             (rt, InputRole.CATALYST)), (switched,))

    cdna_f = g.input("cdna_forward", "10x cDNA-amplification forward primer",
                     MaterialKind.OLIGO, inferred=True,
                     source=SOURCE + " (commercial primer sequence undisclosed)")
    cdna_r = g.input("cdna_reverse", "10x cDNA-amplification reverse primer",
                     MaterialKind.OLIGO, inferred=True,
                     source=SOURCE + " (commercial primer sequence undisclosed)")
    cdna_polymerase = g.input("cdna_polymerase", "cDNA-amplification polymerase",
                              MaterialKind.CATALYST)
    amplified = g.product(
        "amplified", "amplified barcoded host and bacterial cDNA",
        MaterialKind.MOLECULAR, state=barcoded_cdna_scene(), inferred=True,
        source=SOURCE + " (commercial cDNA-amplification primer geometry undisclosed)")
    g.react("amplify_cdna", "Amplify the shared barcoded cDNA pool", StepKind.PCR,
            ((switched, InputRole.TEMPLATE), (cdna_f, InputRole.STRUCTURAL),
             (cdna_r, InputRole.STRUCTURAL),
             (cdna_polymerase, InputRole.CATALYST)), (amplified,),
            geometry=PCRGeometry.unavailable(
                ENRICHED_AMPLICON,
                reason="The protocol uses the commercial 10x cDNA-amplification mix without disclosing one universal primer pair."))
    host = g.product("host_branch", "reserved standard 5-prime GEX aliquot",
                     MaterialKind.MOLECULAR, state=GEX_LIBRARY)
    bacterial = g.product("bacterial_branch", "bacterial 16S enrichment aliquot",
                          MaterialKind.MOLECULAR, state=ENRICHED_AMPLICON)
    g.react("split", "Split host GEX and bacterial-enrichment aliquots",
            StepKind.SPLIT, ((amplified, InputRole.TEMPLATE),),
            (host, bacterial), changes_molecule=False)

    enrich_f = g.input("enrich_forward", "published 16S Enrich Forward primer",
                       MaterialKind.OLIGO)
    enrich_r = g.input("enrich_reverse", "published nested 1061R reverse primer",
                       MaterialKind.OLIGO)
    enrich_polymerase = g.input("enrich_polymerase", "16S enrichment polymerase",
                                MaterialKind.CATALYST)
    first = g.product("first_enrichment", "first-round barcode-bearing 16S amplicons",
                      MaterialKind.MOLECULAR, state=ENRICHED_AMPLICON)
    g.react("enrich_first", "Run the first 35-cycle targeted 16S enrichment",
            StepKind.PCR,
            ((bacterial, InputRole.TEMPLATE), (enrich_f, InputRole.STRUCTURAL),
             (enrich_r, InputRole.STRUCTURAL),
             (enrich_polymerase, InputRole.CATALYST)), (first,),
            geometry=_enrichment_geometry())
    size_rule = g.input("size_rule", "BluePippin 955-1215-bp selection",
                        MaterialKind.SELECTION)
    selected = g.product("selected", "955-1215-bp complete 16S amplicons",
                         MaterialKind.MOLECULAR, state=ENRICHED_AMPLICON)
    off_size = g.product("off_size", "off-size 16S amplification products",
                         MaterialKind.WASTE)
    g.react("size_select", "Retain complete 955-1215-bp 16S products",
            StepKind.SELECTION,
            ((first, InputRole.TEMPLATE), (size_rule, InputRole.SELECTION)),
            (selected, off_size), changes_molecule=False)
    g.discard(off_size, "falls outside the published complete-amplicon size window")
    second = g.product("second_enrichment", "twice-enriched 16S amplicons",
                       MaterialKind.MOLECULAR, state=ENRICHED_AMPLICON)
    g.react("enrich_second", "Run the second 20-cycle targeted 16S enrichment",
            StepKind.PCR,
            ((selected, InputRole.TEMPLATE), (enrich_f, InputRole.STRUCTURAL),
             (enrich_r, InputRole.STRUCTURAL),
             (enrich_polymerase, InputRole.CATALYST)), (second,),
            geometry=_enrichment_geometry())

    reviewed = [
        "associated_cells", "gem_condition", "partition", "poly_dt", "primer_16s",
        "rt", "reverse_transcribe", "bead_oligo", "template_switch",
        "cdna_forward", "cdna_reverse", "cdna_polymerase", "amplify_cdna",
        "split", "enrich_forward", "enrich_reverse", "enrich_polymerase",
        "enrich_first", "size_rule", "size_select", "enrich_second",
    ]
    finish_pcr_branch(
        g, "host", host, GEX_LIBRARY, PRIMERS, reviewed,
        label="INVADEseq host GEX", read_lengths=READ_LENGTHS[
            "Final host gene-expression library"], inferred_primers=True,
        limitation="Dual Index TT Set A is commercial; its lot-specific library-completion primer pair is not printed in the protocol.")
    # Bacterial branch: Steps 111-119 index with Dual Index TT Set A (15 cycles), then
    # repeat the 955-1,215-bp BluePippin program (Steps 92-102) before MiSeq PE300.
    # Written locally because the shared branch helper has no post-PCR selection.
    tt_source = SOURCE + " (endpoint primer sequence not declared in this model)"
    tt_forward = g.input("bacterial_forward", "INVADEseq bacterial 16S forward PCR primer",
                         MaterialKind.OLIGO, inferred=True, source=tt_source)
    tt_reverse = g.input("bacterial_reverse", "INVADEseq bacterial 16S reverse PCR primer",
                         MaterialKind.OLIGO, inferred=True, source=tt_source)
    tt_polymerase = g.input("bacterial_polymerase",
                            "INVADEseq bacterial 16S PCR polymerase", MaterialKind.CATALYST)
    indexed = g.product("bacterial_indexed", "15-cycle dual-indexed 16S products",
                        MaterialKind.MOLECULAR, state=BACTERIAL_LIBRARY, inferred=True,
                        source=tt_source)
    g.react("bacterial_pcr", "Amplify and index the INVADEseq bacterial 16S branch (15 cycles)",
            StepKind.PCR,
            ((second, InputRole.TEMPLATE), (tt_forward, InputRole.STRUCTURAL),
             (tt_reverse, InputRole.STRUCTURAL), (tt_polymerase, InputRole.CATALYST)),
            (indexed,), geometry=PCRGeometry.unavailable(
                BACTERIAL_LIBRARY,
                reason="Dual Index TT Set A is commercial; its lot-specific library-completion primer pair is not printed in the protocol."))
    final_rule = g.input("final_size_rule",
                         "repeat BluePippin 955-1215-bp program (Step 119)",
                         MaterialKind.SELECTION)
    library = g.product("bacterial_library", "INVADEseq bacterial 16S sequencing library",
                        MaterialKind.LIBRARY, state=BACTERIAL_LIBRARY, inferred=True,
                        source=tt_source)
    final_off = g.product("final_off_size", "off-size indexed 16S products",
                          MaterialKind.WASTE)
    g.react("final_size_select", "Re-select indexed 16S products with the same BluePippin program",
            StepKind.SELECTION,
            ((indexed, InputRole.TEMPLATE), (final_rule, InputRole.SELECTION)),
            (library, final_off), changes_molecule=False)
    g.discard(final_off, "falls outside the repeated BluePippin collection window")
    g.sequence(library, SequencingPlan(
        library=BACTERIAL_LIBRARY, primers=PRIMERS,
        acquired_roles=tuple(p.role for p in PRIMERS),
        read_lengths=READ_LENGTHS["Final bacterial 16S enrichment library"]))
    reviewed.extend(("bacterial_forward", "bacterial_reverse", "bacterial_polymerase",
                     "bacterial_pcr", "bacterial_indexed", "final_size_rule",
                     "final_size_select", "bacterial_library", "final_off_size"))
    g.declare_reviewed_scope(reviewed, SOURCE)
    return g.require_publishable()
