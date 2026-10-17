"""GUIDE-seq chemistry model.

Sources read for every sequence in this module
---------------------------------------------
Tsai et al., Nat Biotechnol 33:187-197 (2015), doi:10.1038/nbt.3117

* Online Methods, "dsODN for GUIDE-Seq" -- both modified tag strands.
* Supplementary file MOESM21 (corrected 25 June 2015): Supplementary Fig. 1
  (library cartoon), Supplementary Table 4 (every common primer and adapter),
  and Supplementary Methods, "Protocol for GUIDE-seq Library Preparation"
  (Y-adapter annealing, end repair/A-tailing, ligation, PCR 1, PCR 2, and the
  MiSeq run: custom primers in cartridge positions 13 and 14, cycles
  151 | 8 | 16 | 151 with the paired-end Nextera protocol).

Every ordered library oligo below is published verbatim. The Index 2 reference is
an explicitly qualified computational post-dark-cycle anchor, not a physical
flow-cell oligo ordering sequence. The written forms keep their
modifications; the base strings used for construction are derived from them and
from the canonical ``illumina`` constants, and the ``_require`` checks below fail
at import if a derivation stops agreeing with the printed oligo.
"""
from pathlib import Path
import re
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from batch_ngs import seg
from chemdraw import (Construct, Scene, Segment, bases_pair, complement, feature,
                      oligo, revcomp)
from endprep import dA_tailed_scene, repair_and_dA_tail
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)


def _bases(written: str) -> str:
    """The bases of an oligo written with its modifications and end labels.

    Keeps only the base letters the sources use, so ``[Phos]``, ``P-``, ``*``,
    ``5'-``/``-3'`` and whitespace all drop out without a per-oligo rule.
    """
    return re.sub(r"[^ACGTNW]", "", written)


def _require(ok: bool, message: str) -> None:
    """Fail at construction time, not in a test, when a derivation breaks."""
    if not ok:
        raise ValueError("GUIDE-seq source transcription: " + message)


# --------------------------------------------------------------- the dsODN tag
# Online Methods, "dsODN for GUIDE-Seq": P = 5' phosphate, * = phosphorothioate.
DSODN_A_WRITTEN = "5'-P-G*T*TTAATTGAGTTGTCATATGTTAATAACGGT*A*T-3'"
DSODN_B_WRITTEN = "5'-P-A*T*ACCGTTATTAACATATGACAACTCAATTAA*A*C-3'"
DSODN_A = _bases(DSODN_A_WRITTEN)
DSODN_B = _bases(DSODN_B_WRITTEN)
_require(len(DSODN_A) == 34, "the dsODN is a 34-bp duplex")
_require(DSODN_B == revcomp(DSODN_A),
         "the two written dsODN strands are not complementary")

# ------------------------------------- Supplementary Table 4: P5-side Y adapter
# "Illumina Y-adapters 1-16 (with Molecular Index tag NNWNNWNN)", A01 shown; the
# sixteen members differ only in their 8-nt sample barcode.
A01_WRITTEN = ("AATGATACGGCGACCACCGAGATCTACACTAGATCGCNNWNNWNN"
               "ACACTCTTTCCCTACACGACGCTCTTCCGATC*T")
COMMON_ADAPTER_WRITTEN = "[Phos]GATCGGAAGAGC*C*A"
SAMPLE_BARCODE_A01 = "TAGATCGC"
MOLECULAR_INDEX = "NNWNNWNN"            # Supplementary Methods: N = A/C/T/G, W = A or T
COMMON_BLOCKING_TAIL = "CA"             # 3'-terminal, on two phosphorothioate linkages

Y_ADAPTER_A01 = il.P5 + SAMPLE_BARCODE_A01 + MOLECULAR_INDEX + il.TRUSEQ_READ1
COMMON_ADAPTER = il.STEM + COMMON_BLOCKING_TAIL
_require(Y_ADAPTER_A01 == _bases(A01_WRITTEN),
         "A01 is not P5 + 8-nt barcode + NNWNNWNN + TruSeq Read 1")
_require(COMMON_ADAPTER == _bases(COMMON_ADAPTER_WRITTEN),
         "the MiSeq common adapter is not the TruSeq stem plus a CA tail")
# The 12-bp duplex: the common strand's 5' stem is the reverse complement of the
# Read 1 arm's last 13 bases minus their terminal T, which is the dT overhang.
_require(il.TRUSEQ_READ1.endswith(il.STEM_COMPLEMENT + "T"),
         "the TruSeq Read 1 arm does not end in the adapter stem plus a 3' dT")
_require(revcomp(il.STEM_COMPLEMENT) == il.STEM
         and complement(il.STEM_COMPLEMENT) == il.STEM[::-1],
         "the adapter stem is not self-consistent")
# The remaining two bases of the common strand have no partner on A01: there is no
# second (P7) arm at all. That is the published "single-tail"/half-functional property.
FORK_BASES = il.TRUSEQ_READ1[-15:-13]
_require(len(FORK_BASES) == len(COMMON_BLOCKING_TAIL)
         and not any(bases_pair(a, b)
                     for a, b in zip(COMMON_BLOCKING_TAIL[::-1], FORK_BASES)),
         "the common adapter's 3' CA tail would pair with A01, so it is not a blocking tail")

# ------------------------------- Supplementary Table 4: P7-side arm and primers
# "Custom Sequencing Primers": Read2 is the arm itself, Index1 its reverse complement.
READ2_SEQ_PRIMER = "GTGACTGGAGTCCTCTCTATGGGCAGTCGGTGAT"
INDEX1_SEQ_PRIMER = "ATCACCGACTGCCCATAGAGAGGACTCCAGTCAC"
P7_ARM = READ2_SEQ_PRIMER                  # the AMP-style Read 2 arm, not TruSeq Read 2
_require(INDEX1_SEQ_PRIMER == revcomp(P7_ARM),
         "the custom Index 1 primer is not the reverse complement of the P7 arm")
_require(P7_ARM != il.TRUSEQ_READ2 and P7_ARM != revcomp(il.TRUSEQ_READ2),
         "the P7 arm is not the stock TruSeq Read 2 arm; it must not be modelled as one")

# "P5 Adapter Primers" and the nested gene-specific primers (labels corrected
# 25 June 2015). Each GSP carries a fixed 5' tail plus a dsODN-strand body.
P5_1 = "AATGATACGGCGACCACCGAGATCTA"
P5_2 = "AATGATACGGCGACCACCGAGATCTACAC"
GSP1_TAIL = "GGATCTCGACGCTCTCCCT"
GSP2_TAIL = "CCTCTCTATGGGCAGTCGGTGAT"
GSP1_A = "GGATCTCGACGCTCTCCCTATACCGTTATTAACATATGACA"
GSP1_B = "GGATCTCGACGCTCTCCCTGTTTAATTGAGTTGTCATATGTTAATAAC"
GSP2_A = "CCTCTCTATGGGCAGTCGGTGATACATATGACAACTCAATTAAAC"
GSP2_B = "CCTCTCTATGGGCAGTCGGTGATTTGAGTTGTCATATGTTAATAACGGTA"
GSP2_A_BODY = GSP2_A[len(GSP2_TAIL):]
_require(P5_2 == il.P5 and P5_1 == il.P5[:len(P5_1)],
         "the P5 adapter primers are not nested prefixes of the canonical P5 sequence")
_require(P7_ARM.endswith(GSP2_TAIL),
         "the GSP2 5' tail is not the 3' end of the P7 arm it installs")
_require(GSP1_A[len(GSP1_TAIL):] == DSODN_B[:22]
         and GSP2_A_BODY == DSODN_B[12:]
         and GSP1_B[len(GSP1_TAIL):] == DSODN_A[:29]
         and GSP2_B[len(GSP2_TAIL):] == DSODN_A[6:33],
         "a GSP body is not a verbatim stretch of the dsODN strand it anneals through")

# "P7 Adapters" P701-P708 complete the arm and add the 8-nt P7 barcode.
P701_WRITTEN = ("CAAGCAGAAGACGGCATACGAGATTCGCCTTA"
                "GTGACTGGAGTCCTCTCTATGGGCAGTCGGTGA")
P701_BARCODE = "TCGCCTTA"
_require(P701_WRITTEN == il.P7 + P701_BARCODE + P7_ARM[:-1],
         "P701 is not P7 + 8-nt barcode + the P7 arm less its terminal base")

# ------------------------------------------------------------------- semantics
SAMPLE_BARCODE = feature(
    "guide_seq_y_adapter_barcode", "sample_index", "whitelist",
    whitelist="Tsai et al. 2015 Supplementary Table 4, Y-adapters A01-A16")
P7_BARCODE = feature(
    "guide_seq_p7_barcode", "sample_index", "whitelist",
    whitelist="Tsai et al. 2015 Supplementary Table 4, P7 adapters P701-P708")
MOLECULAR_INDEX_FEATURE = feature(
    "guide_seq_molecular_index", "umi", "degenerate",
    note="NNWNNWNN in the Y adapter: N is A/C/T/G and W is A or T")

TITLE = "GUIDE-seq — cellular capture of nuclease-cut sites"
NOTES = "01_guide-seq.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1038/nbt.3117">Tsai et al., '
          "<i>Nature Biotechnology</i> (2015)</a>, with the corrected (25 June 2015) "
          "supplementary file: Supplementary Fig. 1, Supplementary Table 4 and the "
          "Supplementary Methods library protocol.")
SUMMARY = ("A protected, phosphorylated blunt dsODN integrates into nuclease-created "
           "breaks in living cells. A half-functional Y adapter carrying a sample "
           "barcode and a degenerate molecular index is ligated to the sheared, "
           "dA-tailed ends, and two rounds of nested tag-specific PCR recover the "
           "genomic DNA immediately beside either dsODN end.")
CAVEAT = ("The two STAT-PCR orientations are parallel libraries built from the two dsODN "
          "strands with the +/- GSP pairs. Everything drawn here is the strand-A "
          "orientation; the mirrored reaction is named explicitly at its step. "
          "The recovered PMC NIHMS647337 supplement (Table 4, p. 32) retains the older swapped GSP1 labels. "
          "The corrected publisher supplementary file (notice p. 1) says its original Table 4 "
          "swapped the Nuclease_off_+_GSP1 and Nuclease_off_-_GSP1 labels and omitted "
          "discovery cycling conditions; this model follows the corrected labels and methods. "
          "<b>Published source conflict:</b> the main text, Stage II, says the random "
          "8-bp molecular barcode is incorporated during amplification; the corrected "
          "Supplementary Methods, Supplementary Fig. 1 and Supplementary Table 4 place "
          "the NNWNNWNN index in the Y adapter before ligation. This model installs "
          "that degenerate index at adapter ligation because the printed adapter "
          "oligos locate it there. Whether Stage II means amplification of already "
          "indexed molecules or describes a different step remains unresolved. "
          "<b>Sequencing-model boundary:</b> the manufacturer’s Indexed Sequencing Overview "
          "15057455 v08 (November 2020), p. 4, specifies grafted P5 extension through seven "
          "chemistry-only cycles before the imaged MiSeq Index 2 read. The displayed canonical "
          "P5 anchor marks that effective read-start boundary; it is not a claimed exact "
          "physical grafted primer. The historical instrument-recipe revision is unassigned. "
          "A01 and P701 are printed examples; the final index windows leave experiment-specific "
          "members unselected. "
          "<b>Published source conflict (instrument):</b> the Online Methods name a "
          "\u201cCovaris S200\u201d shearing instrument, while the Supplementary Methods protocol and "
          "its equipment list name the Covaris S2 (S-2 Focused Ultra-sonicator). The model declares "
          "only the source-stated ~500 bp shearing target and does not select an instrument model; "
          "which designation is the typographic error is unresolved here.")

# -------------------------------------------------------------- final library
# Written 5'->3' from P5. The P7-side arm is the GSP2 tail completed by P7_#, so
# this library is a hybrid: TruSeq on the Read 1 side, an AMP-style arm on the
# Read 2 side, which is why the run needs two custom sequencing primers.
DSODN_IN_LIBRARY = revcomp(GSP2_A_BODY)
_require(DSODN_IN_LIBRARY == DSODN_A[:len(GSP2_A_BODY)],
         "the retained dsODN end is not the strand-A prefix the nested primer leaves")

FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"),
    seg("8-nt sample barcode", "S" * len(SAMPLE_BARCODE_A01), "cbc", placeholder=True,
        feature=SAMPLE_BARCODE, length_bp=len(SAMPLE_BARCODE_A01)),
    seg("8-nt molecular index", MOLECULAR_INDEX, "umi", placeholder=True,
        feature=MOLECULAR_INDEX_FEATURE, length_bp=len(MOLECULAR_INDEX)),
    seg("Read 1 arm", il.TRUSEQ_READ1[:-1], "r1"),
    seg("adapter 3' dT : fragment dA", il.TRUSEQ_READ1[-1]),
    seg("adjacent genomic DNA", "X" * 40, placeholder=True),
    seg("integrated dsODN end", DSODN_IN_LIBRARY, "cbc"),
    seg("Read 2 / Index 1 arm", revcomp(P7_ARM), "r2"),
    seg("P7 barcode reverse complement", "I" * len(P701_BARCODE), "cbc", placeholder=True,
        feature=P7_BARCODE, length_bp=len(P701_BARCODE)),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="GUIDE-seq STAT-PCR library")
FINAL_CAPTION = ("One ligated adapter end (P5 side) and one PCR-installed arm (P7 side): "
                 "the dA/dT junction exists only beside the Read 1 arm, and the 8-nt "
                 "degenerate molecular index sits in the adapter ahead of it, not in the "
                 "insert.")

SEQ_PRIMERS = (
    sp.TRUSEQ["R1"],
    sp.custom("Index 1 (i7)", "GUIDE-seq custom Index 1 primer", INDEX1_SEQ_PRIMER,
              "Tsai et al. 2015 Supplementary Table 4; MiSeq cartridge position 13",
              "The exact reverse complement of the GSP2/P7 arm, spiked into the Index "
              "Primer Mix because the stock TruSeq Index 1 primer has no site here."),
    sp.custom("Index 2 (i5)", "effective P5 read-start anchor after seven dark cycles", il.P5,
              'Illumina Indexed Sequencing Overview #15057455 v08 (November 2020), p. 4',
              "MiSeq anneals the template to grafted P5, extends through seven chemistry-only "
              "cycles without imaging, then images the index. The canonical P5 sequence here "
              "is a computational post-extension anchor, not the exact physical grafted oligo. "
              "GUIDE-seq specifies 16 imaged Index 2 cycles; its historical instrument recipe "
              "revision and physical grafted-primer sequence remain unassigned."),
    sp.custom("Read 2", "GUIDE-seq custom Read 2 primer", READ2_SEQ_PRIMER,
              "Tsai et al. 2015 Supplementary Table 4; MiSeq cartridge position 14",
              "The P7 arm itself; Read 2 therefore starts inside the integrated dsODN and "
              "crosses the integration junction into genomic DNA."),
)
# Supplementary Methods: "Sequence with the following number of cycles 151 | 8 | 16 | 151
# with the paired-end Nextera sequencing protocol."
READ_LENGTHS = {"Read 1": 151, "Index 1 (i7)": 8, "Index 2 (i5)": 16, "Read 2": 151}
ACQUIRED_ROLES = tuple(p.role for p in SEQ_PRIMERS)
SEQUENCING_INTRO = ("Read 1 enters the sheared end and reads genomic DNA; Read 2 enters "
                    "the integrated dsODN and crosses the break junction. The 16-cycle "
                    "index-2 read reports the Y-adapter sample barcode and then the "
                    "degenerate molecular index; the 8-cycle index-1 read reports the P7 "
                    "barcode. MiSeq performs seven additional non-imaged chemistry cycles before Index 2; "
                    "those cycles are not added to the source’s 16-cycle reported index read.")


def oligos():
    """Every published oligo, written 5'->3' with its modifications."""
    yield oligo("dsODN strand A — 5' phosphate, phosphorothioate linkages at both ends",
                [Segment("", DSODN_A, "cbc")], mods="P")
    yield oligo("dsODN strand B — 5' phosphate, phosphorothioate linkages at both ends",
                [Segment("", DSODN_B, "cbc")], mods="P")
    yield oligo("Y adapter A01, long strand (A01-A16 differ only in the barcode)", [
        Segment("P5", il.P5, "p5"),
        Segment("8-nt sample barcode", SAMPLE_BARCODE_A01, "cbc"),
        Segment("molecular index", MOLECULAR_INDEX, "umi", placeholder=True),
        Segment("TruSeq Read 1 arm", il.TRUSEQ_READ1[:-1], "r1"),
        Segment("3' dT overhang", il.TRUSEQ_READ1[-1])],
        three=" -3'  (the terminal T sits on a phosphorothioate linkage)")
    yield oligo("MiSeq common adapter, short strand", [
        Segment("adapter stem", il.STEM, "r1"),
        Segment("3' blocking tail", COMMON_BLOCKING_TAIL)],
        mods="[Phos]",
        three=" -3'  (the terminal C and A sit on phosphorothioate linkages)")
    yield oligo("P5_2 nested adapter primer", [Segment("", P5_2, "p5")],
                three=f" -3'  (P5_1 is its first {len(P5_1)} nt)")
    yield oligo("Nuclease_off_+_GSP1", [
        Segment("GSP1 tail", GSP1_TAIL), Segment("dsODN strand B body",
                                                 GSP1_A[len(GSP1_TAIL):], "cbc")])
    yield oligo("Nuclease_off_+_GSP2", [
        Segment("P7 arm tail", GSP2_TAIL, "r2"),
        Segment("dsODN strand B body", GSP2_A_BODY, "cbc")])
    yield oligo("Nuclease_off_-_GSP1", [
        Segment("GSP1 tail", GSP1_TAIL), Segment("dsODN strand A body",
                                                 GSP1_B[len(GSP1_TAIL):], "cbc")])
    yield oligo("Nuclease_off_-_GSP2", [
        Segment("P7 arm tail", GSP2_TAIL, "r2"),
        Segment("dsODN strand A body", GSP2_B[len(GSP2_TAIL):], "cbc")])
    yield oligo("P701 (P701-P708 differ only in the barcode)", [
        Segment("P7", il.P7, "p7"),
        Segment("8-nt P7 barcode", P701_BARCODE, "cbc"),
        Segment("P7 arm", P7_ARM[:-1], "r2")])
    yield oligo("Custom Index 1 sequencing primer", [Segment("", INDEX1_SEQ_PRIMER, "r2")])
    yield oligo("Custom Read 2 sequencing primer", [Segment("", READ2_SEQ_PRIMER, "r2")])


# ------------------------------------------------------------ molecular states
def integration_scene() -> Scene:
    """The dsODN captured blunt at a cellular double-strand break."""
    sc = Scene.duplex([Segment("genome left", "X" * 22, placeholder=True),
                       Segment("dsODN", DSODN_A, "cbc"),
                       Segment("genome right", "X" * 22, placeholder=True)],
                      label="NHEJ-tagged break")
    sc.junction("top", "genome left", "dsODN", "NHEJ capture")
    sc.junction("top", "dsODN", "genome right", "NHEJ capture")
    sc.labels("top")
    return sc


def y_adapter_scene() -> Scene:
    """The annealed half-functional Y adapter, placed by its 12-bp TruSeq stem."""
    top = [seg("P5", il.P5, "p5"),
           seg("8-nt sample barcode", SAMPLE_BARCODE_A01, "cbc", feature=SAMPLE_BARCODE),
           seg("8-nt molecular index", MOLECULAR_INDEX, "umi", placeholder=True,
               feature=MOLECULAR_INDEX_FEATURE, length_bp=len(MOLECULAR_INDEX)),
           seg("Read 1 arm", il.TRUSEQ_READ1[:-13], "r1"),
           seg("adapter stem", il.STEM_COMPLEMENT, "r1"),
           seg("adapter 3' dT overhang", il.TRUSEQ_READ1[-1])]
    bottom = [seg("common adapter stem", il.STEM, "r1"),
              seg("common adapter 3' blocking tail", COMMON_BLOCKING_TAIL)]
    _require(complement(il.STEM_COMPLEMENT) == il.STEM[::-1],
             "the drawn common-adapter stem is not the complement of the A01 stem")
    sc = Scene()
    sc.strand("A01 long strand", top, label="Y adapter A01")
    sc.anneal("MiSeq common adapter", bottom, to="A01 long strand",
              pair=("common adapter stem", "adapter stem"),
              label="MiSeq common adapter", mod5="p",
              unpaired=("common adapter 3' blocking tail",))
    sc.ligation_end("A01 long strand", "adapter 3' dT overhang", kind="cohesive",
                    polarity="3-prime", owner="adapter",
                    text="adapter 3' dT overhang, for the dA-tailed fragment end")
    sc.mark("MiSeq common adapter", "common adapter 3' blocking tail",
            "3'-CA has no partner on A01: there is no second (P7) arm — the "
            "single-tail/half-functional property")
    sc.labels("A01 long strand")
    return sc


def sheared_fragment() -> Construct:
    return Construct([seg("genomic DNA", "X" * 16, placeholder=True),
                      seg("integrated dsODN", DSODN_A, "cbc"),
                      seg("genomic DNA beyond the tag", "X" * 16, placeholder=True)],
                     name="sheared dsODN-tagged genomic fragment")


def dA_tailed_scene_for_fragment() -> Scene:
    return dA_tailed_scene(repair_and_dA_tail(sheared_fragment()),
                           label="end-repaired, dA-tailed fragment")


def adapter_geometry() -> AdapterLigationGeometry:
    """Both ligating ends, as the source prints them."""
    return AdapterLigationGeometry.known(
        "cohesive", y_adapter_scene(), dA_tailed_scene_for_fragment(),
        statement=("Supplementary Table 4 prints both adapter strands: A01 is P5 + an "
                   "8-nt sample barcode + NNWNNWNN + the TruSeq Read 1 arm ending in a "
                   "3' dT, and the MiSeq common adapter is the 12-nt stem plus a "
                   "non-complementary 3' CA block. The Supplementary Methods end-repair "
                   "and A-tail the sheared fragment before ligation, so the 3' dT "
                   "overhang pairs with the fragment's 3' dA."))


def ligation_product_scene() -> Scene:
    """The sealed single-tail capture molecule, with the T:A junction marked."""
    top = [seg("P5", il.P5, "p5"),
           seg("8-nt sample barcode", SAMPLE_BARCODE_A01, "cbc", feature=SAMPLE_BARCODE),
           seg("8-nt molecular index", MOLECULAR_INDEX, "umi", placeholder=True,
               feature=MOLECULAR_INDEX_FEATURE, length_bp=len(MOLECULAR_INDEX)),
           seg("Read 1 arm", il.TRUSEQ_READ1[:-1], "r1"),
           seg("adapter 3' dT : fragment dA", il.TRUSEQ_READ1[-1]),
           seg("genomic DNA", "X" * 16, placeholder=True),
           seg("integrated dsODN", DSODN_A, "cbc"),
           seg("far end", "[genomic DNA and the second ligated adapter]",
               placeholder=True)]
    sc = Scene.duplex(top, label="single-tail capture molecule", unpaired=("far end'",))
    sc.junction("top", "adapter 3' dT : fragment dA", "genomic DNA",
                "adapter-fragment ligation; the fragment's dA is the facing base")
    sc.labels("top")
    return sc


# ------------------------------------------------------------- protocol graph
def _spri(graph, prefix, label, template, beads_label, reason):
    """One source-specified AMPure XP SPRI cleanup with its removed fraction."""
    beads = graph.input(prefix + "_beads", beads_label, MaterialKind.SELECTION)
    kept = graph.product(prefix + "_kept", "SPRI-retained DNA", MaterialKind.MOLECULAR,
                         state=template.state)
    removed = graph.product(prefix + "_removed", "SPRI-removed fraction",
                            MaterialKind.WASTE)
    graph.react(prefix, label, StepKind.CLEANUP,
                ((template, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (kept, removed))
    graph.discard(removed, reason)
    return kept


def protocol_graph():
    graph = ProtocolGraph(TITLE, authority=SOURCE)
    cells = graph.input("cells", "living cells with genomic targets",
                        MaterialKind.BIOLOGICAL, primary=True)
    guide = graph.input("guide_rna", "targeting guide RNA", MaterialKind.OLIGO)
    nuclease = graph.input("programmable_nuclease", "programmable nuclease",
                           MaterialKind.CATALYST)
    cut = graph.product("cut_genome", "cells containing nuclease-created DSBs",
                        MaterialKind.MOLECULAR)
    graph.react("cut", "Create on- and off-target genomic DSBs", StepKind.GENOME_EDITING,
                ((cells, InputRole.TEMPLATE), (guide, InputRole.STRUCTURAL),
                 (nuclease, InputRole.CATALYST)), (cut,))
    dsodn = graph.input("dsodn", "protected phosphorylated blunt 34-bp dsODN tag",
                        MaterialKind.ADAPTER, state=(DSODN_A, DSODN_B))
    nhej = graph.input("nhej", "cellular NHEJ machinery", MaterialKind.CATALYST)
    tagged = graph.product("tagged_genome", "dsODN-tagged genomic break junctions",
                           MaterialKind.MOLECULAR, state=integration_scene())
    graph.react("integrate", "Capture the dsODN at cellular DSBs", StepKind.INTEGRATION,
                ((cut, InputRole.TEMPLATE), (dsodn, InputRole.STRUCTURAL),
                 (nhej, InputRole.CATALYST)), (tagged,))
    dnadvance = graph.input("gdna_beads",
                            "Agencourt DNAdvance solid-phase reversible immobilization beads",
                            MaterialKind.SELECTION)
    gdna = graph.product("genomic_dna", "isolated dsODN-tagged genomic DNA",
                         MaterialKind.MOLECULAR, state=tagged.state)
    cell_debris = graph.product("gdna_removed", "cell lysate fraction not bound to the beads",
                                MaterialKind.WASTE)
    graph.react("isolate_gdna", "Isolate genomic DNA on SPRI beads", StepKind.EXTRACTION,
                ((tagged, InputRole.TEMPLATE), (dnadvance, InputRole.SELECTION)),
                (gdna, cell_debris))
    graph.discard(cell_debris, "is separated from the bead-bound genomic DNA during isolation")
    shear = graph.input("shear_condition", "Covaris shearing to ~500 bp",
                        MaterialKind.CONDITION)
    fragments = graph.product("fragments", "sheared dsODN-tagged genomic fragments",
                              MaterialKind.MOLECULAR, state=sheared_fragment())
    graph.react("shear", "Shear 400 ng of genomic DNA", StepKind.FRAGMENTATION,
                ((gdna, InputRole.TEMPLATE), (shear, InputRole.CONDITION)), (fragments,))
    cleaned = _spri(graph, "shear_cleanup", "Clean the sheared DNA (1X SPRI, elute 15 µL)",
                    fragments, "AMPure XP SPRI beads, 1X ratio",
                    "is removed with the 1X SPRI supernatant after shearing")
    repair_mix = graph.input("end_repair_mix", "end-repair enzyme mix",
                             MaterialKind.CATALYST)
    repaired = graph.product("repaired", "blunt, 5'-phosphorylated fragments",
                             MaterialKind.MOLECULAR)
    graph.react("end_repair", "End-repair the sheared fragments", StepKind.END_REPAIR,
                ((cleaned, InputRole.TEMPLATE), (repair_mix, InputRole.CATALYST)),
                (repaired,))
    taq = graph.input("taq_polymerase", "non-hot-start Taq polymerase (72 C A-tailing)",
                      MaterialKind.CATALYST)
    tailed = graph.product("da_tailed", "end-prepared, 3'-dA-tailed fragments",
                           MaterialKind.MOLECULAR, state=dA_tailed_scene_for_fragment())
    graph.react("a_tail", "Add a single 3' dA to each fragment end", StepKind.TAILING,
                ((repaired, InputRole.TEMPLATE), (taq, InputRole.CATALYST)), (tailed,))
    a01 = graph.input("y_adapter_long_strand",
                      "Y-adapter long strand A01 (P5 + sample barcode + NNWNNWNN + "
                      "TruSeq Read 1 arm + 3' dT)",
                      MaterialKind.ADAPTER, state=Y_ADAPTER_A01)
    common = graph.input("miseq_common_adapter",
                         "5'-phosphorylated MiSeq common adapter with a 3' CA block",
                         MaterialKind.ADAPTER, state=COMMON_ADAPTER)
    y_adapter = graph.product("y_adapter", "annealed half-functional Y adapter",
                              MaterialKind.ADAPTER, state=y_adapter_scene())
    graph.react("anneal_adapter", "Anneal A01 to the MiSeq common adapter",
                StepKind.HYBRIDIZATION,
                ((a01, InputRole.TEMPLATE), (common, InputRole.STRUCTURAL)), (y_adapter,))
    ligase = graph.input("adapter_ligase", "T4 DNA ligase", MaterialKind.CATALYST)
    adapted = graph.product("adapted", "single-tail-adapted genomic fragments",
                            MaterialKind.MOLECULAR, state=ligation_product_scene())
    graph.react("ligate_adapter", "Ligate the half-functional Y adapter to both dA-tailed "
                "ends", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (y_adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (adapted,),
                geometry=adapter_geometry())
    ligated_clean = _spri(graph, "ligation_cleanup",
                          "Clean the ligation (0.9X SPRI, elute 12 µL)", adapted,
                          "AMPure XP SPRI beads, 0.9X ratio",
                          "removes unligated adapter and short products with the 0.9X SPRI supernatant")
    flank_a = graph.product("flank_a", "aliquot for the dsODN strand-A (+) reactions",
                            MaterialKind.MOLECULAR)
    flank_b = graph.product("flank_b", "aliquot for the dsODN strand-B (-) reactions",
                            MaterialKind.MOLECULAR)
    graph.react("split_orientations", "Split into the two junction orientations",
                StepKind.SPLIT, ((ligated_clean, InputRole.TEMPLATE),), (flank_a, flank_b),
                changes_molecule=False)
    p5_1 = graph.input("p5_1", "P5_1 adapter primer", MaterialKind.OLIGO, state=P5_1)
    p5_2 = graph.input("p5_2", "P5_2 nested adapter primer", MaterialKind.OLIGO, state=P5_2)
    gsp1_a = graph.input("gsp1_a", "Nuclease_off_+_GSP1", MaterialKind.OLIGO, state=GSP1_A)
    gsp1_b = graph.input("gsp1_b", "Nuclease_off_-_GSP1", MaterialKind.OLIGO, state=GSP1_B)
    gsp2_a = graph.input("gsp2_a", "Nuclease_off_+_GSP2", MaterialKind.OLIGO, state=GSP2_A)
    gsp2_b = graph.input("gsp2_b", "Nuclease_off_-_GSP2", MaterialKind.OLIGO, state=GSP2_B)
    p7_primer = graph.input("p7_index_primer", "indexed P7 adapter primer (P701-P708)",
                            MaterialKind.OLIGO, state=P701_WRITTEN)
    polymerase = graph.input("pcr_polymerase", "Platinum Taq polymerase",
                             MaterialKind.CATALYST)
    pcr1_reason = ("This declaration does not yet assemble a strand-resolved PCR 1 template "
                   "and product with both half-functional adapter ends. The source prints the "
                   "primers; variable genomic flank length is not proof that their geometry "
                   "cannot be modeled or that source information is absent.")
    pcr2_reason = ("This declaration does not yet assemble the sequential products of GSP2, "
                   "P5_2 and P7_# in the three-primer reaction. The ordered primers are "
                   "source-positive; the unassigned state is a model limitation.")
    amp1_a = graph.product("amp1_a", "strand-A PCR 1 product", MaterialKind.MOLECULAR)
    amp1_b = graph.product("amp1_b", "strand-B PCR 1 product", MaterialKind.MOLECULAR)
    graph.react("stat_pcr1_a", "STAT-PCR 1 from dsODN strand A (GSP1 + P5_1)",
                StepKind.PCR,
                ((flank_a, InputRole.TEMPLATE), (gsp1_a, InputRole.STRUCTURAL),
                 (p5_1, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (amp1_a,), geometry=PCRGeometry.unavailable(FINAL_LIBRARY,
                                                            reason=pcr1_reason))
    graph.react("stat_pcr1_b", "STAT-PCR 1 from dsODN strand B (mirrored GSP1 + P5_1)",
                StepKind.PCR,
                ((flank_b, InputRole.TEMPLATE), (gsp1_b, InputRole.STRUCTURAL),
                 (p5_1, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (amp1_b,), geometry=PCRGeometry.unavailable(FINAL_LIBRARY,
                                                            reason=pcr1_reason))
    clean1_a = _spri(graph, "pcr1_cleanup_a",
                     "Clean strand-A PCR 1 (1.2X SPRI, elute 15 µL; beads carried into PCR 2)",
                     amp1_a, "AMPure XP SPRI beads, 1.2X ratio",
                     "is removed with the 1.2X SPRI supernatant after PCR 1")
    clean1_b = _spri(graph, "pcr1_cleanup_b",
                     "Clean strand-B PCR 1 (1.2X SPRI, elute 15 µL; beads carried into PCR 2)",
                     amp1_b, "AMPure XP SPRI beads, 1.2X ratio",
                     "is removed with the 1.2X SPRI supernatant after PCR 1")
    amp2_a = graph.product("amp2_a", "strand-A nested library", MaterialKind.MOLECULAR,
                           state=FINAL_LIBRARY)
    amp2_b = graph.product("amp2_b", "strand-B nested library", MaterialKind.MOLECULAR)
    graph.react("stat_pcr2_a",
                "Nested STAT-PCR 2 from dsODN strand A (GSP2 + P5_2 + P7_#)", StepKind.PCR,
                ((clean1_a, InputRole.TEMPLATE), (gsp2_a, InputRole.STRUCTURAL),
                 (p5_2, InputRole.STRUCTURAL), (p7_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)),
                (amp2_a,), geometry=PCRGeometry.unavailable(FINAL_LIBRARY,
                                                            reason=pcr2_reason))
    graph.react("stat_pcr2_b",
                "Nested STAT-PCR 2 from dsODN strand B (mirrored GSP2 + P5_2 + P7_#)",
                StepKind.PCR,
                ((clean1_b, InputRole.TEMPLATE), (gsp2_b, InputRole.STRUCTURAL),
                 (p5_2, InputRole.STRUCTURAL), (p7_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)),
                (amp2_b,), geometry=PCRGeometry.unavailable(FINAL_LIBRARY,
                                                            reason=pcr2_reason))
    clean2_a = _spri(graph, "pcr2_cleanup_a",
                     "Clean the strand-A nested library (0.7X SPRI, elute 30 µL)",
                     amp2_a, "AMPure XP SPRI beads, 0.7X ratio",
                     "is removed with the 0.7X SPRI supernatant after PCR 2")
    clean2_b = _spri(graph, "pcr2_cleanup_b",
                     "Clean the strand-B nested library (0.7X SPRI, elute 30 µL)",
                     amp2_b, "AMPure XP SPRI beads, 0.7X ratio",
                     "is removed with the 0.7X SPRI supernatant after PCR 2")
    library = graph.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY)
    graph.react("pool", "Quantify by qPCR, normalize to 1.2e10 molecules and pool both "
                "captured break orientations", StepKind.POOL,
                ((clean2_a, InputRole.TEMPLATE), (clean2_b, InputRole.TEMPLATE)), (library,),
                changes_molecule=False)
    graph.sequence(library, SequencingPlan(library=FINAL_LIBRARY,
                                           primers=tuple(SEQ_PRIMERS),
                                           acquired_roles=ACQUIRED_ROLES,
                                           read_lengths=READ_LENGTHS))
    graph.declare_reviewed_scope((
        "cells", "guide_rna", "programmable_nuclease", "cut", "dsodn", "nhej",
        "integrate", "shear_condition", "shear", "end_repair_mix", "end_repair",
        "taq_polymerase", "a_tail", "y_adapter_long_strand", "miseq_common_adapter",
        "anneal_adapter", "adapter_ligase", "ligate_adapter", "split_orientations",
        "p5_1", "p5_2", "gsp1_a", "gsp1_b", "gsp2_a", "gsp2_b", "p7_index_primer",
        "pcr_polymerase", "stat_pcr1_a", "stat_pcr1_b", "stat_pcr2_a", "stat_pcr2_b",
        "pool", "library", "gdna_beads", "isolate_gdna",
        "shear_cleanup", "shear_cleanup_beads", "ligation_cleanup", "ligation_cleanup_beads",
        "pcr1_cleanup_a", "pcr1_cleanup_a_beads", "pcr1_cleanup_b", "pcr1_cleanup_b_beads",
        "pcr2_cleanup_a", "pcr2_cleanup_a_beads", "pcr2_cleanup_b", "pcr2_cleanup_b_beads"),
        SOURCE)
    return graph.require_publishable()
