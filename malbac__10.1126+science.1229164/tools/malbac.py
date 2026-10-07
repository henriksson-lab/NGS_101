"""MALBAC constructs stated in Zong et al. Science 2012.

The accessible paper specifies the 27-nt common handle and N8 random foot, but the exact
handle bases are confined to an unavailable supplement. The schematic therefore models
those regions honestly as length-preserving placeholders rather than promoting a
secondary-source sequence to a published fact.
"""
from __future__ import annotations

from chemdraw import Construct, Segment, complement_segments

HANDLE_NT = 27
RANDOM_NT = 8
PREAMPLIFICATION_CYCLES = 5
MELT_C = 94
ANNEAL_C = 0
EXTEND_C = 65
LOOP_C = 58


def seg(name: str, seq: str, **kw) -> Segment:
    return Segment(name=name, top=seq, **kw)


def handle(name: str = "common handle") -> Segment:
    return seg(name, "H" * HANDLE_NT, tag="tso", placeholder=True)


def random_foot(name: str = "N8") -> Segment:
    return seg(name, "N" * RANDOM_NT, placeholder=True)


def random_primer() -> Construct:
    return Construct([handle(), random_foot()], name="MALBAC random primer")


def pcr_primer() -> Construct:
    return Construct([handle()], name="MALBAC PCR primer")


def semi_amplicon(insert_nt: int = 28) -> Construct:
    return Construct([handle(), random_foot(),
                      seg("genomic copy", "X" * insert_nt, placeholder=True)],
                     name="MALBAC semi-amplicon")


def full_amplicon(insert_nt: int = 28) -> Construct:
    left = [handle("left handle"), random_foot("left N8"),
            seg("genomic copy", "X" * insert_nt, placeholder=True)]
    return Construct([*left, *complement_segments([handle("right handle")])],
                     name="MALBAC full amplicon")

