"""Shared Tang-lineage whole-transcript cDNA amplification chemistry.

The 2009 supplementary oligo table is the authority for these sequences.  Later
protocols such as scTrio-seq inherit this amplification before attaching a
different sequencing-platform library.
"""
from __future__ import annotations

from chemdraw import Construct, Segment, revcomp


UP1_HANDLE = "ATATGGATCCGGCGCGCCGTCGAC"
UP2_HANDLE = "ATATCTCGAGGGCGCGCCGGATCC"
DT_LEN = 24
UP1 = UP1_HANDLE + "T" * DT_LEN
UP2 = UP2_HANDLE + "T" * DT_LEN

# The second amplification uses the same bases with a 5'-amine modification.
AUP1 = UP1
AUP2 = UP2


def first_strand() -> Construct:
    return Construct([
        Segment("UP1", UP1),
        Segment("first-strand cDNA", "X" * 24, placeholder=True),
    ], name="UP1-primed first-strand cDNA")


def amplified_cdna() -> Construct:
    """Representative ds cDNA after UP1/UP2 amplification."""
    return Construct([
        Segment("UP2 handle", UP2_HANDLE, "r2"),
        Segment("poly(T)", "T" * DT_LEN, "r2"),
        Segment("cDNA", "XXXXXX...XXXXXX", placeholder=True),
        Segment("poly(A)", "A" * DT_LEN, "r1"),
        Segment("UP1 handle'", revcomp(UP1_HANDLE), "r1"),
    ], name="Tang-lineage amplified cDNA")
