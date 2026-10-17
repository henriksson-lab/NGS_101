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
from circular import CircularizationGeometry, circularize_ssdna, rolling_circle
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, StepKind)


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
          "therefore does not assert a fictitious perfectly paired duplex. "
          "The article text (Experimental design, probe design) describes a 23 bp PCR "
          "handle at the probe 5′ end and a 17 bp HindIII-containing extension handle at "
          "the 3′ end, whereas every Supplementary Table 1 probe is ordered as a common "
          "23-nt 5′ flank containing AAGCTT and a common 18-nt 3′ flank, "
          "TGTATGGACCGGGTCATC, that is the exact complement of the 18-nt InDrop reverse "
          "primer. In the final mRNA-complementary orientation the positions agree but "
          "the stated lengths are swapped. The model uses the table sequences: an 18-nt "
          "PCR handle and 23-nt extension handle. The table's variable regions between "
          "these flanks are 55–60 nt, not the 50 bp binding length stated in the text; the "
          "drawing keeps a schematic 50-nt placeholder and does not resolve that "
          "difference.")

# Supplementary Table 1: all 29,765 ordered probes share exactly this 23-nt 5' flank
# and the 18-nt 3' flank below (the reverse complement of the InDrop primer). The
# variable middle is 55-60 nt in the table; it is drawn as a schematic 50-nt
# placeholder (the text's stated binding length) rather than as one fixed pool length,
# so no segment declares TARGET_LENGTH as a known molecular length.
EXTENSION_HANDLE = "AAGCTTTGCTCGAAGAAACTATG"
PCR_HANDLE = "TGTATGGACCGGGTCATC"
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
if revcomp(PCR_HANDLE) != INDROP_REVERSE:
    raise ValueError("the InDrop reverse primer must be the exact PCR-handle complement")
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
            placeholder=True, feature=PROBE_ID),
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


def _phosphorylated_rows() -> list:
    sc = Scene()
    sc.strand("strand", list(ORDERED), label="5′-phosphorylated protoprobe", mod5="p")
    return [*sc.rows(), *annotation_rows(ORDERED)]


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


def extension_template_segments() -> list[Segment]:
    """Supplementary Table 2 'Extension primer (UMI+poly-A)', 5'->3', /3Phos/ blocked."""
    return [
        seg("poly(T)30 template", "T" * 30),
        seg("V+N11 probe-UMI template", "V" + "N" * 11, "umi", placeholder=True,
            feature=PROBE_UMI, length_bp=12),
        seg("extension-handle template", EXTENSION_TEMPLATE_HANDLE, "r2"),
    ]


if "".join(s.top for s in extension_template_segments()) != EXTENSION_TEMPLATE:
    raise ValueError("extension-template segments must spell the printed primer")


def extension_scene() -> Scene:
    sc = Scene()
    sc.strand("protoprobe", list(PROPER), label="RCA-3 protoprobe")
    sc.strand("blocked extension template", extension_template_segments(),
              label="blocked extension template", mod3="p")
    sc.note("blocked extension template",
            "source tables differ by one internal T; strands shown without claimed base alignment")
    sc.arrow("protoprobe", "polymerase copies the 12-base UMI and poly(T)30 template")
    return sc


PRIMER_TABLE = "Supplementary Table 2 (List of Primers)"


FINAL_PROBE = Construct([
    *list(PROPER),
    seg("probe UMI", "U" * 12, "umi", placeholder=True, feature=PROBE_UMI),
    seg("poly(A)30", "A" * 30),
], name="finished ProBac-seq hybridization probe")


def hybridization_scene() -> Scene:
    sc = Scene()
    sc.strand("target transcript", [
        seg("mRNA upstream", "N" * 12, placeholder=True),
        seg("target transcript site", "X" * TARGET_LENGTH, "cbc", placeholder=True),
        seg("mRNA downstream", "N" * 12, placeholder=True),
    ], label="fixed bacterial mRNA")
    sc.anneal("finished probe", list(FINAL_PROBE), to="target transcript",
              pair=("gene-specific probe identity complement", "target transcript site"),
              unpaired=("PCR handle complement", "extension handle / HindIII site complement",
                        "probe UMI", "poly(A)30"), label="hybridized DNA probe")
    return sc


# The ProBac-seq Materials name the 10x Next GEM Single Cell 3' v3.1 kit (PN-1000128,
# v3.1 reagent PN 2000164, Next GEM Chip G) and print the PCR-2 "10X forward" primer
# that matches the 5' end of the gel-bead primer, but print no gel-bead oligo.  The
# 16-nt barcode, 12-nt UMI and T30 below are transcribed from the vendor guide for
# that kit, CG000204 Rev D (fetched into the ProBac bundle by a catalogue
# extra_sources row): its kit list names PN-1000128 and Gel Beads PN 2000164, and its
# Appendix "Oligonucleotide Sequences", Protocol Step 1.5, prints
# CTACACGACGCTCTTCCGATCT-N16-N12-T30VN.  ProBac-seq names no guide revision.
BEAD_GUIDE = ("10x CG000204 Rev D (Chromium Next GEM Single Cell 3′ Reagent Kits v3.1), "
              "Appendix Oligonucleotide Sequences, Gel Bead Primers")
BEAD_NOTE = ("The ProBac-seq Materials name the 10x v3.1 kit (PN-1000128; v3.1 reagent PN "
             "2000164) and print only the 10X forward primer (partial Read 1); the gel-bead "
             "oligo's 16-nt cell barcode, 12-nt UMI and poly(dT)30 are transcribed from that "
             "kit's user guide, " + BEAD_GUIDE + ". ProBac-seq names no guide revision. The "
             "guide's 3′ VN anchor is not drawn as a separate segment: it pairs within the "
             "first two copied probe-UMI positions.")


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
    *bead_oligo()[1:],
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
        oligo("blocked UMI + poly(A) extension template", extension_template_segments(),
              three="/3Phos/-3'"),
        oligo("10x forward primer", [seg("partial Read 1", PARTIAL_R1, "r1")]),
        oligo("InDrop reverse primer", [seg("PCR handle primer", INDROP_REVERSE, "r2")]),
        oligo("PCR-3 forward primer", [seg("P5 + Read 1", PCR3_FORWARD, "p5")]),
        oligo("PCR-3 reverse primer R1", [seg("P7 + i7 + Read 2 + probe handle",
                                                        PCR3_REVERSE, "p7")]),
    ]


def workflow() -> Workflow:
    w = Workflow(MolecularState("Crude ordered Twist protoprobe pool",
                                tuple(_single_strand(ORDERED))))
    w.react("T4 polynucleotide kinase 5′-phosphorylation (37 °C 30 min)",
            _phosphorylated_rows(), name="5′-phosphorylated protoprobe pool",
            note="Step 1: the 5′ phosphate is required for splinted T4 ligation.")
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
    w.react("HindIII-release the final protoprobe pool and gel-purify",
            _single_strand(PROPER), name="final linear protoprobes",
            note="RCA-2 (nicking) primer directs HindIII; 15% TBE–urea gel excision "
                 "(~100–120 nt band) removes primers and short digestion products.")
    w.react("Anneal the 3-prime-blocked extension template",
            extension_scene().rows(), name="extension-template duplex")
    w.react("Klenow-copy the probe UMI and poly(T)30 template, then gel-purify",
            _single_strand(FINAL_PROBE), name="UMI- and poly(A)-tailed probes",
            note="Klenow fragment, room temperature 30 min; a second 15% TBE–urea gel "
                 "separates the unextended blocked primer.")
    w.react("Hybridize probes to fixed bacterial transcripts",
            hybridization_scene().rows(), name="transcript-bound probes",
            note="Stringent washing removes unbound and mismatched probes.")
    w.react("Encapsulate with 10x beads and six-cycle GEM PCR",
            gem_scene().rows(), name="GEM barcode-transfer reaction",
            note=BEAD_NOTE)
    w.react("Recover and PCR-2-amplify cell-barcoded probes",
            Scene.duplex(list(PCR1_PRODUCT), label="PCR-1/PCR-2 product").rows(),
            name="pooled cell-barcoded probe library",
            note="The source-listed partial Read 1 and InDrop primers amplify this product.")
    w.react("PCR-3 add P5, P7 and the sample index",
            Scene.duplex(list(FINAL_LIBRARY), label="indexed library").rows(),
            name="final indexed ProBac-seq library")
    return w


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1038/s41596-024-01002-1"
    graph = ProtocolGraph("ProBac-seq", authority=source)
    crude = graph.input("ordered_protoprobes", "crude ordered Twist protoprobe pool",
                        MaterialKind.MOLECULAR, primary=True, state=ORDERED)
    pnk = graph.input("t4_pnk", "T4 polynucleotide kinase", MaterialKind.CATALYST)
    current = graph.product("phosphorylated_protoprobes", "5′-phosphorylated protoprobe pool",
                            MaterialKind.MOLECULAR, state=_phosphorylated_rows())
    graph.react("phosphorylation", "5′-phosphorylate the oligo pool (Step 1)",
                StepKind.PHOSPHORYLATION,
                ((crude, InputRole.TEMPLATE), (pnk, InputRole.CATALYST)), (current,))
    gel = graph.input("tbe_urea_gel", "15% TBE–urea gel band excision", MaterialKind.SELECTION)
    ligase = graph.input("t4_ligase", "T4 DNA ligase", MaterialKind.CATALYST)
    phi29 = graph.input("phi29", "phi29 DNA polymerase", MaterialKind.CATALYST)
    hindiii = graph.input("hindiii", "HindIII", MaterialKind.CATALYST)
    circles = (CIRCLE1, CIRCLE2, CIRCLE3)
    rcas = (RCA1, RCA2, RCA3)
    primers = (RCA1_PRIMER, RCA2_PRIMER, RCA3_PRIMER)
    for number, (circle, rca, primer_seq) in enumerate(zip(circles, rcas, primers), 1):
        primer = graph.input(f"rca_{number}_primer", f"RCA-{number} splint-primer",
                             MaterialKind.OLIGO, source=PRIMER_TABLE,
                             state=Construct([seg(f"RCA-{number} primer", primer_seq)]))
        closed = graph.product(f"circle_{number}", f"round-{number} splint-closed probe circle",
                               MaterialKind.MOLECULAR, state=circle_rows(circle.linear))
        graph.react(f"circularize_{number}",
                    f"Splint and circularize the round-{number} protoprobe",
                    StepKind.CIRCULARIZATION,
                    ((current, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
                     (ligase, InputRole.CATALYST)), (closed,),
                    geometry=CircularizationGeometry.known(
                        "splinted single-strand ligation",
                        statement=f"The exact RCA-{number} primer spans the future end-to-end ligation bond."))
        concatemer = graph.product(f"concatemer_{number}", f"RCA-{number} concatemer",
                                   MaterialKind.MOLECULAR,
                                   state=_rca_product(rca, number))
        graph.react(f"rca_{number}", f"Phi29 rolling-circle amplify round {number}",
                    StepKind.EXTENSION,
                    ((closed, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
                     (phi29, InputRole.CATALYST)), (concatemer,))
        released = graph.product(f"released_{number}",
                                 f"HindIII-released round-{number} protoprobes",
                                 MaterialKind.MOLECULAR,
                                 state=(PROPER if number in (1, 3) else ORDERED))
        graph.react(f"hindiii_{number}", f"HindIII-release one round-{number} repeat",
                    StepKind.DIGESTION,
                    ((concatemer, InputRole.TEMPLATE), (hindiii, InputRole.CATALYST)),
                    (released,))
        current = released
    final_protoprobes = graph.product("gel_purified_protoprobes",
                                      "gel-purified ~100–120 nt protoprobes",
                                      MaterialKind.MOLECULAR, state=PROPER)
    short_products = graph.product("short_digestion_products",
                                   "nicking primers and short HindIII products",
                                   MaterialKind.WASTE)
    graph.react("protoprobe_gel", "Gel-purify the final protoprobes", StepKind.SELECTION,
                ((current, InputRole.TEMPLATE), (gel, InputRole.SELECTION)),
                (final_protoprobes, short_products), changes_molecule=False)
    graph.discard(short_products, "Steps 20–24 excise only the protoprobe band from the 15% TBE–urea gel.")
    current = final_protoprobes
    extension_template = graph.input("extension_template",
                                     "3′-blocked UMI and poly(T) extension template",
                                     MaterialKind.OLIGO, source=PRIMER_TABLE,
                                     state=Construct(extension_template_segments(),
                                                     name="blocked extension template"))
    annealed = graph.product("extension_duplex", "protoprobe–extension-template duplex",
                             MaterialKind.MOLECULAR, state=extension_scene())
    graph.react("extension_hybridization", "Anneal the blocked extension template",
                StepKind.HYBRIDIZATION,
                ((current, InputRole.TEMPLATE), (extension_template, InputRole.STRUCTURAL)),
                (annealed,))
    extension_polymerase = graph.input("extension_polymerase",
                                      "DNA polymerase I large (Klenow) fragment",
                                      MaterialKind.CATALYST)
    extended = graph.product("extended_probes", "extension reaction products",
                             MaterialKind.MOLECULAR, state=FINAL_PROBE)
    graph.react("probe_extension", "Copy the probe UMI and poly(T)30 template",
                StepKind.EXTENSION,
                ((annealed, InputRole.TEMPLATE), (extension_polymerase, InputRole.CATALYST)),
                (extended,))
    finished = graph.product("finished_probes", FINAL_PROBE.name, MaterialKind.MOLECULAR,
                             state=FINAL_PROBE)
    unextended = graph.product("unextended_primer", "unextended blocked extension primer",
                               MaterialKind.WASTE)
    graph.react("probe_gel", "Gel-purify the UMI- and poly(A)-tailed probes",
                StepKind.SELECTION,
                ((extended, InputRole.TEMPLATE), (gel, InputRole.SELECTION)),
                (finished, unextended), changes_molecule=False)
    graph.discard(unextended, "A 15% TBE–urea gel separates the short unextended blocked primer from the probes.")
    bacteria = graph.input("bacteria", "fixed bacteria containing target transcripts",
                           MaterialKind.BIOLOGICAL, primary=True)
    bound = graph.product("bound_probes", "transcript-bound ProBac probes",
                          MaterialKind.MOLECULAR, state=hybridization_scene())
    graph.react("probe_hybridization", "Hybridize probes to fixed bacterial transcripts",
                StepKind.HYBRIDIZATION,
                ((bacteria, InputRole.TEMPLATE), (finished, InputRole.STRUCTURAL)), (bound,))
    bead = graph.input("gel_bead", "10x 3′ v3.1 gel-bead oligo", MaterialKind.ADAPTER,
                       source=BEAD_GUIDE, state=Construct(bead_oligo(), name="10x 3′ v3.1 gel-bead oligo"))
    gem = graph.product("gem_reaction", "bacterium-and-bead GEM reaction",
                        MaterialKind.MOLECULAR, state=gem_scene())
    excluded = graph.product("partition_loss", "excluded GEM occupancy products",
                             MaterialKind.WASTE)
    graph.react("partition", "Encapsulate bacteria with 10x beads", StepKind.PARTITION,
                ((bound, InputRole.TEMPLATE), (bead, InputRole.STRUCTURAL)), (gem, excluded),
                changes_molecule=False)
    graph.discard(excluded, "GEMs outside the usable bacterium-and-bead occupancy class are excluded")
    # PCR-1 master mix (Procedure) contains only the InDrop reverse primer; the
    # gel-bead oligo is the other primer.
    gem_primers = graph.input("gem_primers", "InDrop reverse primer (PCR-1 master mix)",
                              MaterialKind.OLIGO, source=PRIMER_TABLE,
                              state=Construct([seg("InDrop reverse primer", INDROP_REVERSE)]))
    polymerase = graph.input("pcr_polymerase", "GEM and library PCR polymerase",
                             MaterialKind.CATALYST)
    pcr1 = graph.product("pcr1_product", PCR1_PRODUCT.name, MaterialKind.MOLECULAR,
                         state=PCR1_PRODUCT)
    graph.react("gem_pcr", "Transfer the 10x v3.1 bead barcode, UMI and poly(dT)30 (CG000204 gel-bead primer) by six-cycle GEM PCR", StepKind.PCR,
                ((gem, InputRole.TEMPLATE), (gem_primers, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (pcr1,),
                geometry=PCRGeometry.unavailable(
                    PCR1_PRODUCT, reason="PCR-1 pairs the InDrop reverse primer with the gel-bead oligo; their landing sites are not drawn as one construct.",
                    symmetric_species=True))
    pcr2_primers = graph.input("pcr2_primers", "PCR-2 primer pair: 10X forward and InDrop reverse",
                               MaterialKind.OLIGO, source=PRIMER_TABLE,
                               state=(PARTIAL_R1, INDROP_REVERSE))
    pooled = graph.product("pcr2_product", "pooled cell-barcoded probe library",
                           MaterialKind.MOLECULAR, state=PCR1_PRODUCT)
    graph.react("pcr2", "Recover and PCR-2-amplify cell-barcoded probes", StepKind.PCR,
                ((pcr1, InputRole.TEMPLATE), (pcr2_primers, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (pooled,),
                geometry=PCRGeometry.unavailable(
                    PCR1_PRODUCT, reason="The source-listed partial Read 1 and InDrop primers are represented as one reagent set.",
                    symmetric_species=True))
    p5 = graph.input("pcr3_forward", "PCR-3 P5 forward primer (10X 1step forward)",
                     MaterialKind.OLIGO, source=PRIMER_TABLE,
                     state=Construct([seg("10X 1step forward primer", PCR3_FORWARD)]))
    p7 = graph.input("pcr3_reverse", "PCR-3 indexed P7 reverse primer (10X PCR3 R1 drawn)",
                     MaterialKind.OLIGO, source=PRIMER_TABLE,
                     state=Construct([seg("10X PCR3 R1", PCR3_REVERSE)]))
    library = graph.product("library", FINAL_LIBRARY.name, MaterialKind.LIBRARY,
                            state=FINAL_LIBRARY)
    graph.react("pcr3", "Add P5, P7 and the sample index", StepKind.PCR,
                ((pooled, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY, reason="The exact final construct is modelled, while PCR landing objects are not duplicated."))
    graph.sequence(library, SequencingPlan(
        library=FINAL_LIBRARY, primers=SEQ_PRIMERS, acquired_roles=ACQUIRED_ROLES,
        read_lengths=READ_LENGTHS, read_uses=READ_USES))
    reviewed = ["ordered_protoprobes", "t4_pnk", "phosphorylation", "tbe_urea_gel",
                "t4_ligase", "phi29", "hindiii"]
    for number in range(1, 4):
        reviewed.extend((f"rca_{number}_primer", f"circularize_{number}",
                         f"rca_{number}", f"hindiii_{number}"))
    reviewed.extend((
        "protoprobe_gel", "extension_template", "extension_hybridization",
        "extension_polymerase", "probe_extension", "probe_gel", "bacteria", "probe_hybridization", "gel_bead", "partition",
        "gem_primers", "pcr_polymerase", "gem_pcr", "pcr2_primers", "pcr2",
        "pcr3_forward", "pcr3_reverse", "pcr3", "library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
