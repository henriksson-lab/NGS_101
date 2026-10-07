"""Sequence model for Tang 2009 single-cell mRNA sequencing.

The defining 2009 method used an unbarcoded SOLiD fragment library.  The 16-plex
P2 adapter belongs to the 2010 Nature Protocols extension and is deliberately not
part of :func:`final_library`.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Segment, revcomp


UP1_HANDLE = "ATATGGATCCGGCGCGCCGTCGAC"
UP2_HANDLE = "ATATCTCGAGGGCGCGCCGGATCC"
DT_LEN = 24
UP1 = UP1_HANDLE + "T" * DT_LEN
UP2 = UP2_HANDLE + "T" * DT_LEN

# Amplification uses the same bases with a 5'-amine block.
AUP1 = UP1
AUP2 = UP2

# Unbarcoded SOLiD adapters from the 2009 supplementary oligo table.
P1_TOP = "CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT"
P1_BOTTOM = "ATCACCGACTGCCCATAGAGAGGAAAGCGGAGGCGTAGTGGTT"
P2_TOP = "AGAGAATGAGGAACCCGGGGCAGTT"
P2_BOTTOM = "CTGCCCCGGGTTCCTCATTCTCT"
LIBRARY_PCR_1 = "CCACTACGCCTCCGCTTTCCTCTCTATG"
LIBRARY_PCR_2 = "CTGCCCCGGGTTCCTCATTCT"

# The later 2010 protocol adds a 19-nt internal adapter and one of sixteen 6-nt
# barcodes to P2.  These constants identify the extension without presenting it
# as part of the defining 2009 construct.
BARCODE2010_INTERNAL = "CGCCTTGGCCGTACAGCAG"
BARCODE2010_LEN = 6
BARCODE2010_COUNT = 16


def amplified_cdna() -> Construct:
    """Representative ds cDNA after the second, 5'-amine-blocked PCR."""
    return Construct([
        Segment("UP2 handle", UP2_HANDLE, "r2"),
        Segment("poly(T)", "T" * DT_LEN, "r2"),
        Segment("cDNA", "XXXXXX...XXXXXX", placeholder=True),
        Segment("poly(A)", "A" * DT_LEN, "r1"),
        Segment("UP1 handle'", revcomp(UP1_HANDLE), "r1"),
    ], name="Tang amplified cDNA")


def final_library() -> Construct:
    """2009 unbarcoded SOLiD library, shown in the orientation P1 -> insert -> P2."""
    return Construct([
        Segment("SOLiD P1", P1_TOP, "p5"),
        Segment("cDNA fragment", "XXXXXXXX...XXXXXXXX", placeholder=True),
        # PCR primer 2 leaves this 23-nt P2-derived end in the amplified library.
        Segment("SOLiD P2", revcomp(P2_BOTTOM), "p7"),
    ], name="Tang 2009 unbarcoded SOLiD library")

