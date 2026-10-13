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
MODEL = BranchModel(20, 20, 20, 20)
PREPARATION = "fixed tissue or cells → target retrieval / permeabilization → accessible RNA"
PREPARATION_CAPTION = "Fixation preserves morphology; pretreatment exposes RNA to the probe set. The exact pretreatment depends on specimen and RNAscope kit."
TARGET_CAPTION = "Each probe has an 18–25-nt target-complementary region, spacer and 14-nt Z tail. Only two probes bound adjacently create the complete 28-nt preamplifier site. A standard set contains 20 pairs across roughly 1 kb."
SEQUENCE_BOUNDARY = ("The defining paper withholds the commercial target-probe, "
                     "preamplifier, amplifier and label-probe sequences. The older core "
                     "patent nevertheless discloses exact GAPD double-Z probes and a "
                     "matching preamplifier example; those public sequences are shown "
                     "below but are not asserted to be the current kit reagents.")
DETECTION = "label probes → fluorescent punctum or enzyme-generated chromogenic dot"
DETECTION_CAPTION = "Microscopy localizes individual RNA molecules in their tissue context; label chemistry depends on the fluorescent or bright-field assay."
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
