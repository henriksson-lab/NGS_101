"""
Canonical Illumina and NEBNext oligonucleotide sequences, shared across protocols.

Every sequence here is verbatim from a primary vendor document; the source is named
inline. Nothing in this file is inferred -- protocol-specific guesses belong in the
protocol's own module, flagged there.

Sources
-------
Illumina   "Illumina Adapter Sequences", Document # 1000000002694, v22, Sept 2025
NEB E7805  NEBNext Ultra II FS DNA Library Prep Kit for Illumina, manual v4.0_7/23
NEB E7780  NEBNext Multiplex Oligos, Dual Index Primers Set 2, manual v3.0_6/24
NEB E7600  NEBNext Multiplex Oligos, Dual Index Primers Set 1, manual v6.0_6/24
NEB E7335  NEBNext Multiplex Oligos, Index Primers Sets 1-4, manual v8.0_6/24

Note: Illumina's own document never uses the strings "P5" or "P7" -- those are the
universal community names for the flow-cell lawn oligos.
"""

from __future__ import annotations

from hairpin import Hairpin

# ----------------------------------------------------------- flow cell and read primers
P5 = "AATGATACGGCGACCACCGAGATCTACAC"                    # 29 nt
P7 = "CAAGCAGAAGACGGCATACGAGAT"                         # 24 nt
P5_RC = "GTGTAGATCTCGGTGGTCGCCGTATCATT"                 # as it appears on the bottom strand
P7_RC = "ATCTCGTATGCCGTCTTCTGCTTG"
TRUSEQ_READ1 = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"      # 33 nt
TRUSEQ_READ2 = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"     # 34 nt
INDEX1_PRIMER = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"     # 33 nt, i7 index read
# i5 index read. Reverse-complement workflow (NovaSeq X, NextSeq 1k/2k, iSeq): a dedicated
# primer, = revcomp(TRUSEQ_READ1) (pinned in checks.run_common). Forward-strand workflow
# (MiSeq, HiSeq 2500, NovaSeq 6000 v1.0): no added primer -- the flow-cell P5 oligo primes
# it (Illumina "Indexed Sequencing Overview Guide", Document # 15057455).
INDEX2_PRIMER_RC = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"  # 33 nt

# What you read through into, when the insert is shorter than the read.
# The leading A in each is the dA tail, NOT part of the adapter oligo.
TRIM_SEEN_IN_READ1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
TRIM_SEEN_IN_READ2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# The 12-bp duplex core shared by both arms of a TruSeq-style adapter.
STEM = "GATCGGAAGAGC"
STEM_COMPLEMENT = "GCTCTTCCGATC"

# -------------------------------------------------------------- NEBNext hairpin adaptor
# NEB #E7601A, supplied at 15 uM. Printed in the oligo-kit manuals (E7780 p.2), NOT in
# the library-prep manual. 65 nt: 12-bp stem, 40-nt loop, dU at position 32, 3'-T
# overhang on a phosphorothioate linkage.
#   5'-/5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTC[dU]ACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3'
NEBNEXT_HAIRPIN = ("GATCGGAAGAGCACACGTCTGAACTCCAGTC"      # 1-31, stem + read-2 arm
                   "U"                                    # 32, the uracil, mid-loop
                   "ACACTCTTTCCCTACACGACGCTCTTCCGATCT")   # 33-65, read-1 arm + 3'-T
NEBNEXT_HAIRPIN_DU_POS = 31                                # 0-based index of the U

NEBNEXT_HAIRPIN_MODEL = Hairpin("NEBNext Adaptor for Illumina", NEBNEXT_HAIRPIN,
                                stem_len=12, overhang_3_len=1)
NEBNEXT_USER_OPENED = NEBNEXT_HAIRPIN_MODEL.user_open_at(NEBNEXT_HAIRPIN_DU_POS)

# After USER (UDG + Endo VIII) opens the loop:
NEBNEXT_ARM_READ2 = NEBNEXT_USER_OPENED.left_arm
NEBNEXT_ARM_READ1 = NEBNEXT_USER_OPENED.right_arm

# The full-length read-2-side arm of a canonical TruSeq adapter, for comparison:
# revcomp(TRUSEQ_READ2) minus its leading dA.
CANONICAL_ARM_READ2 = "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"   # 33 nt
# NEB's "truncated design": their read-2 arm is 2 nt short of that, which is why
# libraries built with it need a minimum of 3 PCR cycles to complete the adapters.
NEBNEXT_READ2_TRUNCATION = 2

# ------------------------------------------------------------- NEBNext index primers
# Structure: P7 + [8-nt i7] + TRUSEQ_READ2  /  P5 + [8-nt i5] + TRUSEQ_READ1
# The i7 index as carried in the oligo is the REVERSE COMPLEMENT of the index read.
NEBNEXT_I5_SET1 = {                                        # E7600, i501-i508
    "i501": "TATAGCCT", "i502": "ATAGAGGC", "i503": "CCTATCCT", "i504": "GGCTCTGA",
    "i505": "AGGCGAAG", "i506": "TAATCTTA", "i507": "CAGGACGT", "i508": "GTACTGAC",
}
NEBNEXT_I7_SET1 = {                                        # E7600, i701-i712
    "i701": "CGAGTAAT", "i702": "TCTCCGGA", "i703": "AATGAGCG", "i704": "GGAATCTC",
    "i705": "TTCTGAAT", "i706": "ACGAATTC", "i707": "AGCTTCAG", "i708": "GCGCATTA",
    "i709": "CATAGCCG", "i710": "TTCGCGGA", "i711": "GCGCGAGA", "i712": "CTATCGCT",
}
NEBNEXT_I5_SET2 = {                                        # E7780, i509-i516
    "i509": "TTGCTTGC", "i510": "GAGAGGTT", "i511": "ACCTGGTT", "i512": "AAGCGGAA",
    "i513": "CGGAACAA", "i514": "GGTAAGCT", "i515": "TGTGGCAT", "i516": "ACTACGGA",
}
NEBNEXT_I7_SET2 = {                                        # E7780, i713-i724
    "i713": "AGGAGGAA", "i714": "AGCAAGCA", "i715": "TCATCACC", "i716": "CGTAGGTT",
    "i717": "TCAGATCC", "i718": "CGTGATCA", "i719": "AGTCGCTT", "i720": "GAACGCTT",
    "i721": "TACGCCTT", "i722": "CTCATCAG", "i723": "TCTTCTGC", "i724": "GCTGGATT",
}

# E7335 and friends: the i5 side is an index-less universal primer.
NEBNEXT_UNIVERSAL_PRIMER = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
# Vendor-neutral name for the same 58-mer: P5 + TruSeq Read 1 sharing their ACAC.
TRUSEQ_P5_FULL = NEBNEXT_UNIVERSAL_PRIMER


def nebnext_i5_primer(index: str) -> str:
    """Full NEBNext i5 primer for an 8-nt index, as ordered (5'->3')."""
    return P5 + index + TRUSEQ_READ1


def nebnext_i7_primer(index: str) -> str:
    """Full NEBNext i7 primer for an 8-nt index, as ordered (5'->3')."""
    return P7 + index + TRUSEQ_READ2


# --------------------------------------------------------------------- NEB thermal data
# E7805 fragmentation time -> insert size. Pick the row, don't guess.
NEBNEXT_FS_SIZING = {           # minutes at 37 C -> (min bp, max bp)
    30: (100, 250),
    20: (150, 350),
    15: (200, 450),
    10: (300, 700),
    5: (500, 1000),
}

# NOT PUBLISHED, recorded so nobody goes looking:
#   NEBNext UMI adaptor sequence -- NEB deleted it from the manual at rev 3.0 (9/25)
#   as "inaccurate" and declared the correct sequence proprietary. Any copy found in an
#   older manual or on a third-party page is unreliable. Do not use one.
