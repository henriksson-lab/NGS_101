"""Source-backed molecular model for ProBac-seq probe capture."""
from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from chemdraw import (Construct, MolecularState, Scene, Segment, Workflow,
                      annotation_rows, circle_rows, complement_segments,
                      feature, oligo, revcomp)
from circular import circularize_ssdna, rolling_circle


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


TITLE = "ProBac-seq — probe-based bacterial single-cell transcript counting"
NOTES = "01_probac-seq.html"
SOURCE = ('<a href="https://doi.org/10.1038/s41596-024-01002-1">'
          'Samanta et al., <i>Nature Protocols</i> (2024)</a>.')
SUMMARY = ("A splint closes each gene-specific protoprobe for three alternating rounds "
           "of rolling-circle amplification. A blocked template then installs a probe UMI "
           "and poly(A), and 10x GEM PCR adds a cell barcode to probes bound in bacteria.")
CAVEAT = ("Supplementary Tables 1 and 2 differ by one T in the extension-handle region: "
          "the ordered probes print AAGCTTTGCT…, while the blocked template prints "
          "AAGCTTGCT…. Both are retained verbatim below; the extension-template panel "
          "therefore does not assert a fictitious perfectly paired duplex.")

# Supplementary Table 1. The middle, target-specific region is represented at the
# protocol's nominal 50 nt rather than pretending that all members of the pool are equal.
EXTENSION_HANDLE = "AAGCTTTGCTCGAAGAAACTATG"
PCR_HANDLE = "GTATGGACCGGGTCATC"
TARGET_LENGTH = 50

# Supplementary Table 2, exactly as written (the extension template has V + 11 N =
# twelve degenerate positions and a 3-prime phosphate).
RCA1_PRIMER = "TTCGAGCAAAGCTTGATGACCCGGTCC"
RCA2_PRIMER = "GGACCGGGTCATCAAGCTTTGCTCGAA"
RCA3_PRIMER = RCA1_PRIMER
EXTENSION_TEMPLATE_HANDLE = "AAGCTTGCTCGAAGAAACTATG"
EXTENSION_TEMPLATE = "T" * 30 + "V" + "N" * 11 + EXTENSION_TEMPLATE_HANDLE
PARTIAL_R1 = "CTACACGACGCTCTTCCGATCT"
INDROP_REVERSE = "GATGACCCGGTCCATACA"
PCR3_FORWARD = il.TRUSEQ_P5_FULL
PCR3_I7 = "AGTAAACC"                 # Supplementary Table 2 primer R1 example.
PCR3_REVERSE = il.P7 + PCR3_I7 + il.TRUSEQ_READ2 + INDROP_REVERSE

CELL = feature("cell_barcode", "cell_barcode", "whitelist",
               whitelist="10x Chromium 3-prime v3 whitelist")
GEM_UMI = feature("gem_umi", "umi", "random",
                  note="12-nt UMI supplied by the 10x gel-bead oligo")
PROBE_UMI = feature("probe_umi", "umi", "random",
                    note="12 degenerate positions copied from the blocked template")
PROBE_ID = feature("probac_probe", "feature_barcode", "whitelist",
                   whitelist="designed ProBac-seq probe pool")
SAMPLE = feature("sample_index_i7", "sample_index", "fixed",
                 note="one of the sixteen Supplementary Table 2 i7 primers")


def ordered_protoprobe() -> Construct:
    """One representative member of the ordered Twist pool, 5-prime to 3-prime."""
    return Construct([
        seg("extension handle / HindIII site", EXTENSION_HANDLE, "r2"),
        seg("gene-specific probe identity", "X" * TARGET_LENGTH, "cbc",
            placeholder=True, feature=PROBE_ID, length_bp=TARGET_LENGTH),
        seg("PCR handle", PCR_HANDLE, "r1"),
    ], name="ordered protoprobe")


ORDERED = ordered_protoprobe()
PROPER = Construct(complement_segments(list(ORDERED), suffix=" complement"),
                   name="mRNA-complementary protoprobe")


def _same_circle(a: str, b: str) -> bool:
    # ``complement_segments`` lower-cases schematic placeholders whereas the
    # sequence-only RCA primitive leaves them upper-case. They denote the same
    # unknown target region and are irrelevant to the circular rotation check.
    a, b = a.upper(), b.upper()
    return len(a) == len(b) and b in a + a


def _rca(template: Construct, primer: str):
    circle = circularize_ssdna(template, five_prime_phosphate=True)
    return circle, rolling_circle(circle, primer, copies=2)


CIRCLE1, RCA1 = _rca(ORDERED, RCA1_PRIMER)
CIRCLE2, RCA2 = _rca(PROPER, RCA2_PRIMER)
CIRCLE3, RCA3 = _rca(ORDERED, RCA3_PRIMER)
if not _same_circle(RCA1.unit, PROPER.top()):
    raise ValueError("RCA-1 plus HindIII must release the complementary protoprobe")
if not _same_circle(RCA2.unit, ORDERED.top()):
    raise ValueError("RCA-2 plus HindIII must restore the ordered orientation")
if not _same_circle(RCA3.unit, PROPER.top()):
    raise ValueError("RCA-3 plus HindIII must release the hybridizing orientation")


def _single_strand(con: Construct, label: str = "") -> list:
    sc = Scene()
    sc.strand("strand", list(con), label=label or con.name)
    return [*sc.rows(), *annotation_rows(con)]


def _splint_scene(template: Construct, primer: str, round_no: int) -> Scene:
    """Draw the exact primer across the last-to-first circularization junction."""
    site = revcomp(primer)
    boundary = template.top()[-len(PCR_HANDLE):] + template.top()[:len(EXTENSION_HANDLE)]
    hits = [i for i in range(len(boundary) - len(site) + 1)
            if boundary.startswith(site, i)]
    if len(hits) != 1:
        raise ValueError(f"RCA-{round_no} primer does not uniquely bridge the probe ends")
    sc = Scene()
    sc.strand("probe ends", [seg("joined-end landing site", site, "r2")],
              label=f"round {round_no} probe ends brought together")
    sc.anneal("splint / RCA primer", [seg(f"RCA-{round_no} primer", primer, "r1")],
              to="probe ends", pair=(f"RCA-{round_no} primer", "joined-end landing site"),
              label=f"RCA-{round_no} splint and primer")
    sc.mark("probe ends", "joined-end landing site", "primer spans the future ligation bond")
    return sc


def _rca_product(product, round_no: int) -> Construct:
    unit = product.unit
    return Construct([
        seg(f"RCA-{round_no} complementary unit 1", unit, "r2", placeholder=True,
            length_bp=len(unit)),
        seg(f"RCA-{round_no} complementary unit 2", unit, "r2", placeholder=True,
            length_bp=len(unit)),
    ], name=f"RCA-{round_no} concatemer (two units drawn)")


def extension_scene() -> Scene:
    sc = Scene()
    sc.strand("protoprobe", list(PROPER), label="RCA-3 protoprobe")
    template = [
        seg("poly(T)30 template", "T" * 30),
        seg("V+N11 probe-UMI template", "V" + "N" * 11, "umi", placeholder=True,
            feature=PROBE_UMI, length_bp=12),
        seg("extension-handle template", EXTENSION_TEMPLATE_HANDLE, "r2"),
    ]
    sc.strand("blocked extension template", template, label="blocked extension template",
              mod3="p")
    sc.note("blocked extension template",
            "source tables differ by one internal T; strands shown without claimed base alignment")
    sc.arrow("protoprobe", "polymerase copies the 12-base UMI and poly(T)30 template")
    return sc


FINAL_PROBE = Construct([
    *list(PROPER),
    seg("probe UMI", "U" * 12, "umi", placeholder=True, feature=PROBE_UMI),
    seg("poly(A)30", "A" * 30),
], name="finished ProBac-seq hybridization probe")


def hybridization_scene() -> Scene:
    sc = Scene()
    sc.strand("target transcript", [
        seg("mRNA upstream", "N" * 12, placeholder=True),
        seg("target transcript site", "X" * TARGET_LENGTH, "cbc", placeholder=True,
            length_bp=TARGET_LENGTH),
        seg("mRNA downstream", "N" * 12, placeholder=True),
    ], label="fixed bacterial mRNA")
    sc.anneal("finished probe", list(FINAL_PROBE), to="target transcript",
              pair=("gene-specific probe identity complement", "target transcript site"),
              unpaired=("PCR handle complement", "extension handle / HindIII site complement",
                        "probe UMI", "poly(A)30"), label="hybridized DNA probe")
    return sc


def bead_oligo() -> list[Segment]:
    return [
        seg("partial Read 1", PARTIAL_R1, "r1"),
        seg("cell barcode", "B" * 16, "cbc", placeholder=True, feature=CELL),
        seg("10x UMI", "U" * 12, "umi", placeholder=True, feature=GEM_UMI),
        seg("poly(dT)30", "T" * 30),
    ]


def gem_scene() -> Scene:
    sc = Scene()
    sc.strand("probe", list(FINAL_PROBE), label="probe still bound to bacterial mRNA")
    sc.anneal("released gel-bead primer", bead_oligo(), to="probe",
              pair=("poly(dT)30", "poly(A)30"),
              unpaired=("partial Read 1", "cell barcode", "10x UMI"),
              label="one gel-bead oligo")
    sc.arrow("released gel-bead primer", "six-cycle GEM PCR copies and cell-barcodes the probe")
    return sc


PCR1_PRODUCT = Construct([
    *bead_oligo(),
    seg("copied probe UMI", "U" * 12, "umi", placeholder=True, feature=PROBE_UMI),
    *list(ORDERED),
], name="cell-barcoded probe after GEM PCR")


FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"),
    seg("Read 1 arm remainder", il.TRUSEQ_P5_FULL[len(il.P5):], "r1"),
    seg("cell barcode", "B" * 16, "cbc", placeholder=True, feature=CELL),
    seg("10x UMI", "U" * 12, "umi", placeholder=True, feature=GEM_UMI),
    seg("poly(dT)30", "T" * 30),
    seg("probe UMI", "U" * 12, "umi", placeholder=True, feature=PROBE_UMI),
    *list(ORDERED),
    seg("Index 1 / Read 2 arm", revcomp(il.TRUSEQ_READ2), "r2"),
    seg("i7 reverse complement", revcomp(PCR3_I7), "cbc", feature=SAMPLE),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="indexed ProBac-seq library")

SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
_problems = sp.verify(FINAL_LIBRARY, SEQ_PRIMERS,
                      required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
if _problems:
    raise ValueError("invalid ProBac-seq final library: " + "; ".join(_problems))

FINAL_CAPTION = ("One representative probe-counting amplicon. The two UMIs have distinct "
                 "origins: the first is on the 10x gel-bead oligo; the second was copied "
                 "from the blocked extension template before hybridization.")
SEQUENCING_INTRO = ("The published run acquires Index 1 for 8 cycles and Read 1 for 119 "
                    "cycles. The final molecule also contains the standard TruSeq Read 2 "
                    "primer site, although Read 2 is not acquired in that run. Read 1 is "
                    "split computationally and the probe UMI replaces the GEM UMI for "
                    "molecule counting.")
READ_LENGTHS = {"Read 1": 119, "Index 1 (i7)": 8}
ACQUIRED_ROLES = ("Read 1", "Index 1 (i7)")
READ_USES = {"Read 1": "extract", "Index 1 (i7)": "demultiplex"}
SHOW_READ_CONTENT = True


def oligos() -> list[str]:
    return [
        oligo("representative ordered protoprobe architecture", list(ORDERED)),
        oligo("first/third-round RCA splint-primer",
              [seg("RCA-1/RCA-3", RCA1_PRIMER, "r1")]),
        oligo("second-round RCA splint-primer", [seg("RCA-2", RCA2_PRIMER, "r1")]),
        oligo("blocked UMI + poly(A) extension template", [
            seg("poly(T)30", "T" * 30),
            seg("V+N11", "V" + "N" * 11, "umi", placeholder=True, length_bp=12),
            seg("extension handle", EXTENSION_TEMPLATE_HANDLE, "r2")],
              three="/3Phos/-3'"),
        oligo("10x forward primer", [seg("partial Read 1", PARTIAL_R1, "r1")]),
        oligo("InDrop reverse primer", [seg("PCR handle primer", INDROP_REVERSE, "r2")]),
        oligo("PCR-3 forward primer", [seg("P5 + Read 1", PCR3_FORWARD, "p5")]),
        oligo("PCR-3 reverse primer R1", [seg("P7 + i7 + Read 2 + probe handle",
                                                        PCR3_REVERSE, "p7")]),
    ]


def workflow() -> Workflow:
    w = Workflow(MolecularState("Ordered protoprobe pool",
                                tuple(_single_strand(ORDERED))))
    w.react("Anneal the RCA-1 primer across both probe ends",
            _splint_scene(ORDERED, RCA1_PRIMER, 1).rows(),
            name="RCA-1 splinted probe",
            note="The same 27-mer both aligns the ends for ligation and primes phi29.")
    w.react("T4-ligate the phosphorylated protoprobe into a circle",
            circle_rows(CIRCLE1.linear, "RCA-1 splint-directed closure"),
            name="first-round probe circle")
    w.react("Phi29 rolling-circle amplify",
            _single_strand(_rca_product(RCA1, 1)), name="RCA-1 concatemer",
            note="Two complementary repeat units are drawn; the physical concatemer is longer.")
    w.react("Anneal RCA-2 and HindIII-release one complementary unit",
            [*_splint_scene(PROPER, RCA2_PRIMER, 2).rows(), *_single_strand(PROPER)],
            name="mRNA-complementary protoprobe",
            note="HindIII acts at the duplex made by the primer over the AAGCTT site.")
    w.react("Re-circularize and phi29-amplify with RCA-2",
            [*circle_rows(CIRCLE2.linear, "RCA-2 splint-directed closure"),
             *_single_strand(_rca_product(RCA2, 2))], name="RCA-2 concatemer",
            note="This round reverses orientation again.")
    w.react("HindIII-release, re-circularize and amplify with RCA-3",
            [*_splint_scene(ORDERED, RCA3_PRIMER, 3).rows(),
             *circle_rows(CIRCLE3.linear, "RCA-3 splint-directed closure"),
             *_single_strand(_rca_product(RCA3, 3))], name="RCA-3 concatemer",
            note="The third round restores the mRNA-complementary output orientation.")
    w.react("HindIII-release the final protoprobe pool",
            _single_strand(PROPER), name="final linear protoprobes")
    w.react("Anneal the 3-prime-blocked extension template",
            extension_scene().rows(), name="extension-template duplex")
    w.react("Isothermally copy the probe UMI and poly(T)30 template",
            _single_strand(FINAL_PROBE), name="UMI- and poly(A)-tailed probes")
    w.react("Hybridize probes to fixed bacterial transcripts",
            hybridization_scene().rows(), name="transcript-bound probes",
            note="Stringent washing removes unbound and mismatched probes.")
    w.react("Encapsulate with 10x beads and six-cycle GEM PCR",
            gem_scene().rows(), name="GEM barcode-transfer reaction")
    w.react("Recover and PCR-2-amplify cell-barcoded probes",
            Scene.duplex(list(PCR1_PRODUCT), label="PCR-1/PCR-2 product").rows(),
            name="pooled cell-barcoded probe library",
            note="The source-listed partial Read 1 and InDrop primers amplify this product.")
    w.react("PCR-3 add P5, P7 and the sample index",
            Scene.duplex(list(FINAL_LIBRARY), label="indexed library").rows(),
            name="final indexed ProBac-seq library")
    return w
