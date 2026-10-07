"""
CRISPR-MIP: sgRNA quantification by molecular inversion probe.

Selinger M, Yakovenko I, Nazir I, Henriksson J. "CRISPR-MIP replaces PCR and reveals GC
and oversampling bias in pooled CRISPR screens." bioRxiv 2024.03.28.587082.
doi:10.1101/2024.03.28.587082

The idea: replace the readout PCR of a pooled CRISPR screen with a padlock capture. The
probe carries a UMI, so each captured genomic molecule is labelled before any
amplification -- which turns the readout from "count reads" into "count molecules", and
removes the PCR-derived GC and sequencing-depth biases the paper reports.
"""

from __future__ import annotations

import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Segment, revcomp
from crispr import SCAFFOLD_V1
from padlock import Padlock

# --------------------------------------------------------------- probe design
# All sequences verbatim from Table S2 of the preprint (supplementary media-2.xlsx).
EXT_ARM = "GTGGAAAGGACGAAACACC"        # 19 nt, published Tm 53.8 C -- ends at the U6 +1
LIG_ARM = "AGCTAGGTCTTGAAAGGAGTGGG"    # 23 nt, published Tm 58.3 C -- past the scaffold

UMI = "N" * 13                                             # 13 nt, published
# The backbone's two sequencing-primer sites are the Illumina TruSeq ones, referenced from
# lib/illumina.py rather than retyped. Table S2 (ref/TableS2_primers_and_probes.xlsx; the
# selftest checks every probe against it) shows the per-probe variable stretch is 8 nt, so
# the constant Read 1 site is TruSeq Read 1 less its 5' ACAC.
READ2_SITE = revcomp(il.TRUSEQ_READ2)                      # 34 nt
READ1_SITE = il.TRUSEQ_READ1[4:]                           # 29 nt
INDEX_LEN = 8

# The nine probes differ ONLY in this 8-nt i7 index (Table S2, MIP_probe_1..9).
PROBE_INDICES = ("ATCACGAC", "ACAGTGGT", "CAGATCCA", "ACAAACGG", "ACCCAGCA",
                 "AACCCCTC", "CCCAACCT", "CACCACAC", "GAAACCCA")

# An example Brunello-style (not G-initiated) spacer, for drawings and tests.
EXAMPLE_SPACER = "ATCGATCGATCGATCGATCG"

CAPTURED_SPAN = 112      # published, and computed from the real vector
UMI_LEN = 13
READ1_CYCLES = 60        # published: reads the captured sgRNA
READ2_CYCLES = 15        # published: enough for the 13-nt UMI

# The amplification primers. These are SEPARATE oligos, and -- the key design point --
# they anneal inside the CAPTURED sequence (the tracrRNA/scaffold), not in the probe
# backbone. A probe that annealed but was never extended therefore has no primer sites.
# Their 5' ends are the Illumina P5/P7 flow-cell sequences.
P5_TRACR_TAIL = "GTGCTTTTTTGAATTCGC"          # Table S2: vector just past the scaffold
P7_TRACR_TAIL = revcomp(SCAFFOLD_V1[32:50])   # Table S2: GTTGATAACGGACTAGCC, on the scaffold
P5_TRACR_FWD = il.P5 + P5_TRACR_TAIL
P7_TRACR_REV = il.P7 + P7_TRACR_TAIL


def backbone(index: str = PROBE_INDICES[0]) -> str:
    """The 84-nt non-complementary middle of the probe."""
    return UMI + READ2_SITE + index + READ1_SITE


PROBE_LEN = len(LIG_ARM) + len(backbone()) + len(EXT_ARM)   # 126
BACKBONE_LEN = len(backbone())                              # 84

# The paper's Methods text says the probe is "134 bp". Table S2's actual sequences are
# 126 nt, and decompose exactly as above. We follow Table S2 and flag the discrepancy.
PROBE_LEN_IN_TEXT = 134

def probe_construct(index: str = PROBE_INDICES[0]) -> Construct:
    """The probe as ordered, 5'->3': lig arm - UMI - Read 2 site - i7 - Read 1 site - ext arm.

    One definition, used by the page diagram, the selftest and the notes (via mdfacts), so
    a length or an index width cannot be restated anywhere and go stale.
    """
    return Construct([
        Segment("ligation arm", LIG_ARM, "r3"),
        Segment("UMI", UMI, "umi", placeholder=True),
        Segment("Read 2 site", READ2_SITE, "t7"),
        Segment("i7", index, None),
        Segment("Read 1 site", READ1_SITE, "s5"),
        Segment("extension arm", EXT_ARM, "r1"),
    ], name="CRISPR-MIP probe")


def probe(index: str = PROBE_INDICES[0]) -> Padlock:
    """One of the nine CRISPR-MIP probes, selected by its i7 index."""
    return Padlock(f"CRISPR-MIP-{index}", ext_arm=EXT_ARM, lig_arm=LIG_ARM,
                   backbone=backbone(index))


# ------------------------------------------------------------------ protocol
# Verbatim conditions from the paper's four steps.
PROTOCOL = (
    ("1. denature & hybridise",
     "1-10 ug gDNA, 0.2-0.002 uM probe, 1.5x Ampligase buffer. 94 C 5 min, then ramp to "
     "60 C at -0.1 C/s, hybridise overnight at 60 C."),
    ("2. extend & ligate",
     "Add Ampligase + 4 U Phusion HF + 0.2 mM dNTPs, pre-heated to 60 C. 60 C, 1 h."),
    ("3. exonuclease",
     "37 C, 45 min with 10 U Exonuclease I (removes un-circularised probe) and 50 U "
     "Exonuclease III (removes genomic DNA). Inactivate 80 C, 20 min."),
    ("4. PCR",
     "Amplify the whole reaction. The P5/P7 primer sites lie WITHIN the captured "
     "sequence, so only probes that were both extended and ligated amplify."),
)

# Why the method is specific, in one place: three independent filters, none of which is a
# size selection or a cleanup.
SELECTION_STEPS = (
    "both arms must anneal to the same molecule, a fixed distance apart",
    "the polymerase must cross the gap and the ligase must close the circle",
    "exonuclease I/III destroy everything still linear -- probe, gDNA, failed captures",
    "the PCR primer sites sit inside the captured region, not the backbone",
)


# ------------------------------------------------------------ sequencing primers
# Sequences come from lib/ via seqprimers; where each lands, and what it reads first, is
# computed from the final library on the page. The library has no i5 index: P5 runs
# straight into captured vector sequence.
SEQ_PRIMERS = [
    replace(sp.TRUSEQ["R1"], note="Its 5' ACAC sits over the last bases of the 8-nt i7, "
            "so how many of those 4 pair varies by probe; the 3' 29 nt pair on all nine."),
    sp.TRUSEQ["I1"],
    sp.mismatching(sp.TRUSEQ["I2"],
                   "its 3' ...GTGT needs the ACAC that starts a full TruSeq Read 1 site; "
                   "the probe has only the 3' 29 nt of that site, preceded by the i7's last "
                   "bases, so the primer's 3' end is mismatched -- and there is no i5 to read"),
    sp.custom("Index 2 (i5)", "P5 flow-cell oligo (forward-strand i5 workflow)", il.P5,
              'Illumina "Indexed Sequencing Overview Guide" #15057455',
              "No i5 in this library: these cycles read constant captured-vector sequence."),
    sp.TRUSEQ["R2"],
]
