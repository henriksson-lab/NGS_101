"""
Segment definitions and constants for Small-seq -- single-cell small-RNA sequencing.

Source
------
Hagemann-Jensen M, Abdullayev I, Sandberg R, Faridani OR.
"Small-seq for single-cell small-RNA sequencing." Nat Protoc 2018;13(10):2407-2424.
doi:10.1038/s41596-018-0049-y          (extracted text: pdf/hagemann.txt, gitignored)
Original method: Faridani OR et al., Nat Biotechnol 2016;34(12):1264-1266.

Why this chemistry is in the repo
--------------------------------
It is the only protocol here that is *not* built on TruSeq DNA adapters and not built on
either template switching or tagmentation. Four things are new:

  1. **TruSeq Small RNA adapters** (RA5 / RA3) rather than the TruSeq DNA adapter, so the
     read-1 landing site is a completely different sequence (see 02_sequencing_primer.md).
  2. **Sequential ligation** -- 3' adapter first, then a digestion step that destroys the
     unligated 3' adapter, then the 5' adapter. No polymerase is involved in adapter
     attachment at all.
  3. **A UMI in the 5' adapter**, not in a TSO and not on a probe. It is the first thing
     sequenced, before any of the insert.
  4. **A masking oligonucleotide** that occludes the 3' end of 5.8S rRNA so that the 3'
     adapter cannot be ligated to it. Depletion by hybridisation and *inaction*, rather
     than by pulldown.

Evidence marking (repo convention)
----------------------------------
  [green]  verbatim from the paper
  [yellow] derived or inferred here
  [red]    not published -- not in the paper's main text at all
  [check]  computed and asserted in tools/selftest.py

Every real sequence below is [green] from Reagent setup (p. 2414), which states
"All oligonucleotides are listed in the 5' to 3' direction". Line-wrap spaces in the
extracted text have been closed up; nothing else was changed.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Segment, feature, revcomp
from illumina import P5, P7, P7_RC
import seqprimers as sp

UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "whitelist",
                     whitelist="small-seq-srx-192-index-set")
FEATURES = {"UMI": UMI_FEATURE, "i7": I7_FEATURE}

# ======================================================================= oligos
# TODO(lib): move RA5 / RA3 / RTP (TruSeq Small RNA adapters) to illumina.py.
# All [green] verbatim, Reagent setup. Modifications are carried separately as _MODS
# strings so that the base sequence stays computable (tm(), revcomp(), len()).

# --- 5' adapter: RNA, 5'-amino-blocked, free 3'-OH, UMI + CA at the 3' end -----------
# RA5 (NH2-rGrUrUrCrArGrArGrUrUrCrUrArCrArGrUrCrCrGrArCrGrArUrCrHrHrHrHrHrHrHrHrCrA)
RA5_HANDLE = "GTTCAGAGTTCTACAGTCCGACGATC"      # 26 nt, written as DNA for arithmetic
RA5_UMI_LEN = 8                                 # [green] "a stretch of eight 'H'"
RA5_UMI_BASE = "H"                              # [green] rH = rA, rU or rC
RA5_UMI_DIVERSITY = 3 ** RA5_UMI_LEN            # [green] "6,561 different UMIs"
RA5_LINKER = "CA"                               # [green] rCrA, separates UMI from insert
RA5_MODS = "NH2 (aminolink C6) at the 5' end; 3'-OH; every base a ribonucleotide"

# --- 3' adapter: DNA, 5'-adenylated, 3'-blocked ---------------------------------------
# RA3 (rAppTGGAATTCTCGGGTGCCAAGG-ddC)
RA3 = "TGGAATTCTCGGGTGCCAAGG"                   # 21 nt
RA3_MODS = "5'-rApp (pre-adenylated); 3'-ddC"

# --- RT primer, which is also the reverse primer of PCR 1 -----------------------------
# RTP (biotin-CCTTGGCACCCGAGAATTCCrA)
RTP = "CCTTGGCACCCGAGAATTCCA"                   # 21 nt
# Kept verbatim although it equals revcomp(RA3): the paper prints the two independently, and
# that equality (asserted in selftest) is what the digestion step depends on.
RTP_MODS = "5'-biotin (blocks lambda exonuclease); 3'-terminal rA"

# --- forward primer of both PCRs ------------------------------------------------------
# RP1 (AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA)   <- as printed, 50 nt
# Built from its parts so that P5 (lib/illumina.py) and RA5 are not typed twice; the
# selftest asserts the result equals the printed string in the comment above.
RP1_ANNEAL_NT = 21          # [check] its 3' 21 nt == RA5_HANDLE[:21]; the 5' 29 are P5
RP1 = P5 + RA5_HANDLE[:RP1_ANNEAL_NT]           # 50 nt

# --- 5.8S rRNA masking oligonucleotide ------------------------------------------------
# [green] ATCGGCAAGCGACGCTCAGA CAGGCGTAGCCCCGGG AGGAACCCGGGGC CGCAAGTGCGTTCGAAGTGTCGATGAT-biotin
MASK_58S = ("ATCGGCAAGCGACGCTCAGA"
            "CAGGCGTAGCCCCGGG"
            "AGGAACCCGGGGC"
            "CGCAAGTGCGTTCGAAGTGTCGATGAT")      # 76 nt
MASK_58S_MODS = "3'-biotin, to block extension"
# [yellow] what it must be annealing to. The paper says only that it "targets the 3' end
# of the rRNA molecule" and "binds to the 3' end of the highly abundant 5.8S rRNA"; it
# never prints the rRNA. Checking this base-by-base needs a reference sequence, which is
# data and therefore does not live in this repo (README, "Data policy"). Recorded so the
# next person can BLAST it rather than retype it.
MASK_58S_TARGET = revcomp(MASK_58S)
RRNA_58S_NT = 156           # [yellow] textbook length of mature human 5.8S; NOT in the paper

# --- index primers: SEQUENCE NOT PUBLISHED in the main text ---------------------------
# [red] "SRX DNA index primers are modified from standard Illumina TruSeq small RNA index
# primers to have barcodes of 8 bp for 192 samples. For sequence information, see
# Supplementary Table 1." -- and Supplementary Table 1 is not in the extracted text.
# Do not invent it. What IS published:
SRX_SEQUENCE = None                             # [red] deliberately absent
SRX_INDEX_NT = 8                                # [green] "barcodes of 8 bp"
SRX_N_INDICES = 192                             # [green] "for 192 samples" (two 96 plates)
SRX_CONC_UM = 50                                # [green] dissolved to 50 uM, 1 ul per well

# The right-hand arm it must add to the library, beyond RA3. Open question as a switch
# (README, "Open questions are switches, not prose"): flip SRX_LINKER_NT and every size
# on the page re-derives.
#
#   SRX 3' end must anneal where RTP anneals, i.e. over RA3 -- that length is already in
#   the construct and adds nothing. What is *added* is P7 + index + whatever links them
#   to the RTP-like 3' end:
SRX_P7_NT = len(P7)                             # 24 [check] = len(illumina.P7)
SRX_LINKER_NT = 12                              # [yellow] see derivation below
SRX_RIGHT_ARM_NT = SRX_P7_NT + SRX_INDEX_NT + SRX_LINKER_NT      # 44

# Derivation of SRX_LINKER_NT, from two independent published sizes:
#   known left arm   = P5 + RA5 + UMI + CA + RA3 = 86 bp   [check]
#   (a) miRNA window [green] "aiming to cut targets in the size range of 145-155 bp",
#       a typical miRNA being 22 nt -> right arm = 145-86-22 .. 155-86-22 = 37..47
#   (b) 5.8S peak    [green] "a peak appears at the 270- to 290-bp range ... if the
#       rRNA-masking oligonucleotide is omitted"; 5.8S ~156 nt [yellow]
#                           -> right arm = 270-86-156 .. 290-86-156 = 28..48
#   intersection 37..47. A TruSeq-small-RNA-shaped arm (P7 + 8-bp index + 12-nt link)
#   is 44, which sits inside both windows -- asserted in selftest as a bracketing, not
#   as a sequence.

# ================================================================= protocol constants
LYSIS_TEMP_C, LYSIS_MIN = 72, 20                # [green] Step 9
LIG3_TEMP_C, LIG3_HOURS = 30, 6                 # [green] Step 12, then 4 C for 10 h
LIG5_TEMP_C, LIG5_MIN = 37, 60                  # [green] Step 18
RT_TEMP_C, RT_MIN = 42, 60                      # [green] Step 21, SuperScript II
PCR1_ANNEAL_C, PCR1_CYCLES = 60, 13             # [green] Step 24, cycles 2-14
PCR2_ANNEAL_C, PCR2_CYCLES = 67, 13             # [green] Step 31, cycles 2-14
PCR1_PRIMERS = ("RP1", "RTP (leftover from the digestion step)")   # [green]
PCR2_PRIMERS = ("RP1", "SRX index primer")                         # [green]

READ_LEN = 51                                   # [green] Step 38, "a single read of 51 bp"
TRIM_5P = RA5_UMI_LEN + len(RA5_LINKER)         # 10 nt removed before mapping [green]
SMALLRNA_MAX_NT = 40                            # [green] maxRlen, "40 bp or shorter"
PRECURSOR_MIN_NT = 41                           # [green] minRlen, reads of maximum length
SMALLRNA_MIN_NT = 18                            # [green] cutadapt --minimum-length 18

PIPPIN_WINDOW_BP = (130, 160)                   # [green] Step 35A(i) start/end
GEL_CUT_BP = (120, 200)                         # [green] Step 35B(iv)
SIZE_SELECT_TARGET_BP = (145, 155)              # [green] "aiming to cut targets in"
PIPPIN_RESULT_BP = (140, 180)                   # [green] Anticipated results, Fig. 3d
RRNA_PEAK_BP = (270, 290)                       # [green] Fig. 3b, mask omitted
LIBRARY_PROFILE_BP = (100, 300)                 # [green] Anticipated results, Fig. 3a

TYPICAL_MIRNA_NT = 22                           # [yellow] not stated; used for sizing only

# Enzymes, in order. [green] Reagent tables of Steps 10, 13, 16, 19, 22.
ENZYMES = [
    ("3' ligation", "T4 RNA ligase 2, truncated KQ",
     "no ATP, so only the pre-adenylated adapter can be ligated -- small RNAs cannot be "
     "ligated to each other"),
    ("digestion", "5' deadenylase",
     "strips the 5'-rApp from unligated RA3, leaving a 5'-phosphate"),
    ("digestion", "lambda exonuclease",
     "digests the 5'-phosphorylated strand of a duplex: RTP anneals to free RA3 and the "
     "exonuclease eats the RA3, while RTP's 5'-biotin protects it"),
    ("5' ligation", "T4 RNA ligase 1 + ATP",
     "needs a 5'-phosphate on the acceptor, which capped mRNA does not have"),
    ("RT", "SuperScript II, in Taq buffer",
     "Taq buffer rather than RT buffer because excess MgCl2 is carried over"),
    ("PCR", "Phusion Hot Start II", "high fidelity, both rounds"),
]


# ======================================================================= constructs
def _seg(name, top, tag=None, **kw):
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


def ra5_oligo() -> str:
    """RA5 as ordered, written as DNA (H = rA/rU/rC placeholder)."""
    return RA5_HANDLE + RA5_UMI_BASE * RA5_UMI_LEN + RA5_LINKER


def insert_placeholder(n: int) -> str:
    """A stand-in insert of n characters -- never complemented as if it were bases."""
    if n <= 8:
        return "X" * n
    return "X" * (n - 6) + "...XXX"


def ligated_rna(insert_nt: int = TYPICAL_MIRNA_NT) -> Construct:
    """Small RNA after both ligations: RA5 + UMI + CA + insert + RA3. Single-stranded."""
    return Construct([
        _seg("RA5", RA5_HANDLE, "r1"),
        _seg("UMI", RA5_UMI_BASE * RA5_UMI_LEN, "umi", placeholder=True),
        _seg("CA", RA5_LINKER, "tso"),
        _seg("small RNA", insert_placeholder(insert_nt), None, placeholder=True),
        _seg("RA3", RA3, "r2"),
    ], name="ligated small RNA")


def library(insert_nt: int = TYPICAL_MIRNA_NT) -> Construct:
    """The final, indexed library.

    Left of RA3 every base is published. Right of RA3 nothing is: SRX lives in
    Supplementary Table 1, so the linker and the index are drawn as placeholders and
    marked inferred. P7 has to be there for the library to cluster at all.
    """
    return Construct([
        _seg("P5", P5, "p5"),
        _seg("RA5", RA5_HANDLE, "r1"),
        _seg("UMI", RA5_UMI_BASE * RA5_UMI_LEN, "umi", placeholder=True),
        _seg("CA", RA5_LINKER, "tso"),
        _seg("small RNA", insert_placeholder(insert_nt), None, placeholder=True),
        _seg("RA3", RA3, "r2"),
        _seg("SRX link", "L" * SRX_LINKER_NT, "r3", placeholder=True, inferred=True,
             note="length derived from the published library sizes; sequence not published"),
        _seg("i7", "A" * SRX_INDEX_NT, "cbc", placeholder=True, inferred=True,
             note="8 bp, 192 of them; sequences in Supplementary Table 1"),
        _seg("P7", P7_RC, "p7", inferred=True,
             note="required for clustering; SRX itself is not published"),
    ], name="Small-seq library")


def library_bp(insert_nt: int) -> int:
    """Library length for a given insert, from the segment arithmetic alone."""
    return (len(P5) + len(RA5_HANDLE) + RA5_UMI_LEN + len(RA5_LINKER)
            + insert_nt + len(RA3) + SRX_RIGHT_ARM_NT)


LEFT_ARM_NT = len(P5) + len(RA5_HANDLE) + RA5_UMI_LEN + len(RA5_LINKER) + len(RA3)   # 86
ADAPTER_DIMER_BP = library_bp(0)                 # [yellow] 130 bp -- see 01_small-seq.md


# ============================================================== sequencing primers
# By reference: sequences live in lib/ (TruSeq DNA) or above (the Small RNA oligos). The
# section on the page computes where each lands on library() and what it reads first.
#
# TODO(lib): move to illumina.py -- the TruSeq Small RNA block (RA5, RA3, RTP, RP1, RPI
# and the Small RNA Sequencing / Index primers), verbatim from the Illumina adapter
# document. Until then the read-1 entry below is the site the CONSTRUCT requires (RA5, 26
# nt, ending where the UMI begins), not a vendor primer sequence.
READ1_SITE = RA5_HANDLE

SEQ_PRIMERS = [
    sp.custom("Read 1", "TruSeq Small RNA read-1 site (RA5) -- required by the construct",
              READ1_SITE, "derived: [green] \"Sequencing starts from the UMI\"",
              "The paper names no sequencing primer. Any usable one is a suffix of P5 + RA5 "
              "ending on RA5's last base; the vendor Small RNA read-1 primer is not yet in "
              "lib/ to check against."),
    sp.custom("Read 1", "RP1 re-used as read-1 primer (counter-example)", RP1,
              "Reagent setup, p. 2414",
              "Matches, but stops 5 nt short of RA5's 3' end, so read 1 would start CGATC "
              "instead of on the UMI -- contradicting the paper. Not what was used."),
    sp.mismatching(sp.TRUSEQ["R1"], "no TruSeq DNA adapter in this library; a run loaded "
                                    "only with it reads nothing"),
    sp.mismatching(sp.TRUSEQ["I1"], "the i7 read must prime from RA3 into the SRX linker, "
                                    "which is unpublished (Supplementary Table 1); there is "
                                    "no TruSeq DNA Index 1 site"),
    sp.mismatching(sp.TRUSEQ["I2"], "single-indexed: RP1 puts P5 directly on RA5, so there "
                                    "is no i5 index and no TruSeq i5 primer site"),
    sp.mismatching(sp.TRUSEQ["R2"], "single-read protocol (Step 38, \"a single read of 51 "
                                    "bp\"); no TruSeq DNA read-2 site"),
]
