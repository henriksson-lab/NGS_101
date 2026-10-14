"""Canonical SMART-seq2 oligos and reusable cDNA/PCR geometry.

Protocol pages that explicitly inherit SMART-seq2 should construct the branch here
instead of restating its oligos or calling them unavailable.
"""
from __future__ import annotations

from chemdraw import Construct, Segment, revcomp
import nextera as nx
from pcr import PCRGeometry, PCRPrimer
import rt

DT_LEN = 30
DT_ANCHOR = "VN"
DT_LINKER = "AC"
TSO_LINKER = "ACAT"
OLIGO_DT = rt.oligo_dt(DT_LEN, DT_ANCHOR, rt.SMART_HANDLE + DT_LINKER)
TSO = rt.SMART_HANDLE + TSO_LINKER + rt.TSO_G_TAIL_LNA
ISPCR = rt.SMART_HANDLE


def preamp_template(insert_nt: int = 36) -> Construct:
    """Full-length cDNA after template switching and before symmetric ISPCR."""
    if insert_nt < 1:
        raise ValueError("SMART-seq2 cDNA insert must be non-empty")
    return Construct([
        Segment("5′ SMART handle", rt.SMART_HANDLE, "tso"),
        Segment("TSO linker", TSO_LINKER),
        Segment("TSO G contribution", "GGG", "tso"),
        Segment("transcript cDNA", "X" * insert_nt, placeholder=True),
        Segment("poly(A) copy", "A" * DT_LEN),
        Segment("oligo-dT linker reverse complement", revcomp(DT_LINKER)),
        Segment("3′ SMART handle reverse complement", revcomp(rt.SMART_HANDLE), "tso"),
    ], name="SMART-seq2 full-length cDNA before ISPCR")


def ispcr_geometry(template: Construct | None = None) -> PCRGeometry:
    template = template or preamp_template()
    return PCRGeometry.symmetric(
        template, PCRPrimer("ISPCR", ISPCR),
        statement=("SMART-seq2 places the same 23-nt ISPCR handle at both cDNA ends, "
                   "so one primer species supplies both inward-facing sites."))


def nextera_index_geometry(final_library: Construct) -> PCRGeometry:
    """Invariant N5xx/N7xx landing sites; representative indices do not affect binding."""
    return PCRGeometry.exact(
        final_library,
        PCRPrimer("N5xx", nx.n5xx_primer("A" * 8), anneal_nt=len(nx.S5)),
        PCRPrimer("N7xx", nx.n7xx_primer("A" * 8), anneal_nt=len(nx.S7)),
        statement=("Nextera index primers terminate in s5 and s7; those 3′ regions face "
                   "inward on unlike-ended tagmented cDNA. Representative A8 index "
                   "windows show the invariant landing geometry."))
