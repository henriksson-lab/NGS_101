"""RNAscope double-Z branched-DNA signal-amplification architecture."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from branched_dna import BranchModel, SequenceExample, paired_tail_docking_scene
from chemdraw import Segment, revcomp
from topology import AdjacentProbeDockingTopology

TITLE = "RNAscope"
NOTES = "01_rnascope.html"
SOURCE = ('Defining source: <a href="https://doi.org/10.1016/j.jmoldx.2011.08.002">'
          'Wang et al., <i>Journal of Molecular Diagnostics</i> (2012)</a>; architecture '
          'cross-checked against current Bio-Techne RNAscope documentation. A public '
          'double-Z implementation comes from '
          '<a href="https://patents.google.com/patent/US7709198B2/en">US7709198B2</a>.')
SUMMARY = "Adjacent double-Z probes recognize RNA in fixed cells or tissue and nucleate a branched-DNA signal tree for fluorescent or chromogenic microscopy. The RNA itself is not amplified."
# Wang 2012 Methods (RNAscope Design Strategy): "typically" 20 pairs per ~1 kb, but
# "We chose 10 to 20 pairs"; 20 amplifier sites per preamplifier, 20 label sites per amplifier.
MODEL = BranchModel(10, 20, amplifiers_per_preamp=20, labels_per_amplifier=20)
TREE_MULTIPLICITY_SOURCE = ("Wang et al. 2012 (Methods, RNAscope Design Strategy) state 20 amplifier "
                            "binding sites per preamplifier and 20 label-probe binding sites "
                            "per amplifier.")
PREPARATION = "fixed cells or FFPE section → (FFPE: citrate retrieval) → protease → accessible RNA"
PREPARATION_CAPTION = ("Wang 2012 Methods: cultured cells on slides are fixed in 4% formaldehyde for 60 min and digested with protease (2.5 µg/mL) at 23–25 °C. 5-µm FFPE sections are deparaffinized in xylene, taken through an ethanol series, boiled (100–103 °C) in citrate buffer pH 6 for 15 min (printed as \"10 nmol/L\"; a nanomolar retrieval buffer is implausible and 10 mmol/L is likely intended, but the model does not assert the corrected unit), rinsed and digested with 10 µg/mL protease at 40 °C for 30 min. Current commercial kits may use different pretreatment reagents.")
TARGET_CAPTION = "Each probe has an 18–25-nt target-complementary region, spacer and 14-nt Z tail. Only two probes bound adjacently create the complete 28-nt preamplifier site. The paper typically targets ~1 kb with 20 pairs and used 10–20 pairs; three pairs already gave visible HRP/DAB signal. Target probes hybridize at 40 °C for 3 h in buffer A (6× SSC, 25% formamide, 0.2% lithium dodecyl sulfate, blocking reagents), followed by three room-temperature washes (0.1× SSC, 0.03% LDS)."
SEQUENCE_BOUNDARY = ("The defining paper withholds the commercial target-probe, "
                     "preamplifier, amplifier and label-probe sequences. The older core "
                     "patent nevertheless discloses exact GAPD double-Z probes and a "
                     "matching preamplifier example; those public sequences are shown "
                     "below but are not asserted to be the current kit reagents.")
DETECTION = "label probes (Alexa Fluor or HRP) → fluorescent punctum or DAB dot"
DETECTION_CAPTION = ("Wang 2012 Methods: label probes were conjugated to Alexa Fluor 488, 546, 647 or 750 for epifluorescence, or to HRP for DAB chromogen with hematoxylin counterstain and bright-field microscopy (the design section also names alkaline phosphatase/Fast Red). Spots are counted per cell; UBC and bacterial dapB probe sets serve as positive and negative controls.")
HYBRIDIZATION_SERIES = ("Wang 2012 Methods: sequential 40 °C hybridizations, each followed by three room-temperature washes (0.1× SSC, 0.03% LDS): preamplifier 2 nM in buffer B (20% formamide, 5× SSC, 0.3% LDS, 10% dextran sulfate, blocking reagents) for 30 min; amplifier 2 nM in buffer B for 15 min; label probe 2 nM in buffer C (5× SSC, 0.3% LDS, blocking reagents) for 15 min. Multiplex assays use equimolar reagents of each amplification system.")
READOUT = "RNAscope produces spatial microscopy signal, not a sequencing library. There are therefore no flow-cell adapters, sample indices, sequencing primers or nucleotide reads to draw."

# US7709198B2, Table 2 and SEQ ID NO:1.  Each GAPD oligo is target-binding
# sequence + five-thymidine spacer + one half of the cooperative docking site.
GAPD217_TARGET = "CCAGTGGACTCCACGACGTAC"
GAPD218_TARGET = "TTCTCCATGGTGGTGAAGACG"
SPACER = "TTTTT"
Z_TAIL_1 = "GAAGTTACCGTTTT"
Z_TAIL_2 = "CTGAGTCAAAGCAT"
AMPLIFIER_SITE = "AGGCATAGGACCCGTGTCT"
# The patent prints this terminal 28-mer explicitly. Deriving it from the two probe
# tails makes the cooperative pairing an invariant of the construct rather than a
# duplicated transcription that can silently drift.
PREAMPLIFIER_DOCKING = revcomp(Z_TAIL_1 + Z_TAIL_2)
SEQUENCE_TOPOLOGY = AdjacentProbeDockingTopology(
    "patent-disclosed GAPD double-Z junction",
    revcomp(GAPD217_TARGET),
    revcomp(GAPD218_TARGET),
    SPACER,
    SPACER,
    Z_TAIL_1,
    Z_TAIL_2,
    PREAMPLIFIER_DOCKING,
)


def _probe(name: str, target: str, tail: str) -> tuple[Segment, ...]:
    return (Segment(f"{name} target-binding region", target, "me"),
            Segment(f"{name} T5 spacer", SPACER),
            Segment(f"{name} Z tail", tail, "cbc"))


def preamplifier() -> tuple[Segment, ...]:
    return (Segment("amplifier site 1", AMPLIFIER_SITE, "r2"),
            Segment("T10 spacer", "T" * 10),
            Segment("amplifier site 2", AMPLIFIER_SITE, "r2"),
            Segment("T5 spacer", SPACER),
            Segment("double-Z docking site", PREAMPLIFIER_DOCKING, "r1"))


def sequence_example() -> SequenceExample:
    return SequenceExample(
        "Patent-disclosed double-Z sequence example",
        "US7709198B2 reports these oligos for a human GAPD experiment. They establish "
        "one exact implementation of the cooperative 28-base junction, not the identity "
        "of present-day commercial RNAscope reagents.",
        (("GAPD217 double-Z probe", _probe("GAPD217", GAPD217_TARGET, Z_TAIL_1)),
         ("GAPD218 double-Z probe", _probe("GAPD218", GAPD218_TARGET, Z_TAIL_2)),
         ("example preamplifier (SEQ ID NO:1)", preamplifier())),
        paired_tail_docking_scene(Z_TAIL_1, Z_TAIL_2, PREAMPLIFIER_DOCKING),
        "The docking site is computed and checked as the reverse complement of the two "
        "14-nt tails. The tails reside on separate target probes; the shared line only "
        "makes their cooperative duplex with the preamplifier legible.")
