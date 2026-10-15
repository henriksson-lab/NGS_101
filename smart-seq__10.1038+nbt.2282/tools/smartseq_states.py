"""SMART-handle molecular scenes shared by related SMART-seq pages.

This module owns molecular states and drawing geometry. Each protocol assembles its
own reaction graph and declares its own source evidence.
"""
from __future__ import annotations

import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
sys.path.insert(0, str(HERE.parents[1] / "smart-seq2__10.1038+nmeth.2639" / "tools"))

import nextera as nx
import rt
import smartseq as ss
import smartseq2 as s2
from chemdraw import Construct, Scene, Segment, complement_segments, revcomp

H = rt.SMART_HANDLE

def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)

def ss2_tso_linker(protocol: str) -> list[Segment]:
    """TSO linker segments, 5'->3', for one SMART-handle member.

    FLASH-seq's supplement prints AC outright. SMART-seq2 gets AC plus an inferred AT.
    SMART-seq's own TSO is an undisclosed commercial oligo, so its whole linker is
    inferred.
    """
    if protocol == "FLASH-seq":
        return [seg("AC", s2.DT_LINKER)]
    if protocol == "SMART-seq":
        return [seg("ACAT", s2.TSO_LINKER, inferred=True)]
    return s2.defining_tso_linker()

BODY = 20           # drawn transcript body
POLYA = 30          # drawn poly(A): exactly the oligo-dT's (T)30, so the handle overhangs


def mrna() -> list[Segment]:
    return [seg("body", "X" * BODY, placeholder=True),
            seg("B", "B", placeholder=True),
            seg("polyA", "A" * POLYA)]


def _dt_primer(handle: Segment, link: list[Segment]) -> list[Segment]:
    """Anchored oligo-dT, 5'->3': handle, (linker), T30, V, N."""
    return [handle, *link, seg("dT", "T" * POLYA),
            seg("V", "V", placeholder=True), seg("N", "N", placeholder=True)]


def _cdna(handle: Segment, link: list[Segment]) -> list[Segment]:
    """First strand 5'->3': the oligo-dT, the copied body (the N already paired the last
    body base), then MMLV's untemplated CCC."""
    return [*_dt_primer(handle, link),
            *complement_segments([seg("body", "X" * (BODY - 1), placeholder=True)]),
            seg("CCC", rt.UNTEMPLATED_TAIL, "tso")]


def _rt_scenes(handle: Segment, link: list[Segment], tso: list[Segment]) -> tuple[Scene, ...]:
    """(priming, first strand + CCC, template switch) for one chemistry."""
    s1 = Scene()
    s1.strand("mRNA", mrna())
    s1.anneal("oligo-dT", _dt_primer(handle, link), to="mRNA", pair=("dT", "polyA"))
    s1.mark("oligo-dT", "N", "VN anchor: sits at the poly(A) junction", through="V")
    s1.arrow("oligo-dT", "MMLV")

    s2 = Scene()
    s2.strand("mRNA", mrna())
    s2.anneal("cDNA", _cdna(handle, link), to="mRNA", pair=("dT", "polyA"))
    s2.mark("cDNA", "CCC", "untemplated CCC, overhanging the mRNA 5' end")

    s3 = Scene()
    s3.strand("mRNA", mrna())
    s3.anneal("cDNA", _cdna(handle, link), to="mRNA", pair=("dT", "polyA"))
    s3.anneal("TSO", tso, to="cDNA", pair=("rGrG", "CCC"), above=True)
    s3.mark("cDNA", "CCC", "CCC pairs the TSO's G tail")
    s3.arrow("cDNA", "RT switches template")
    return s1, s2, s3

SS2_HANDLE = seg("H", H, "tso")


def ss2_tso(protocol: str = "SMART-seq2") -> list[Segment]:
    """The drawn TSO for one SMART-handle member, 5'->3'. The G tail is spelt GGG."""
    if protocol == "SMART-seq2":
        return s2.drawn_tso()
    inf = protocol == "SMART-seq"
    return [seg("H", H, "tso", inferred=inf), *ss2_tso_linker(protocol),
            seg("rGrG", "GGG", "tso", inferred=inf)]


SS2_TSO = ss2_tso()                                                    # rGrG+G, as GGG


def ss2_rt(protocol: str = "SMART-seq2") -> tuple[Scene, ...]:
    return _rt_scenes(SS2_HANDLE, [seg("AC", s2.DT_LINKER)], ss2_tso(protocol))

def ss2_preamp_template(protocol: str = "SMART-seq2") -> Construct:
    if protocol == "SMART-seq2":
        return s2.drawn_preamp_template(BODY, POLYA)
    inf = protocol == "SMART-seq"
    return Construct([
        seg("H5", H, "tso", inferred=inf),
        *[seg("TSO linker" + ("" if i == 0 else f" {s.name}"), s.top, s.tag,
              inferred=s.inferred or inf)
          for i, s in enumerate(ss2_tso_linker(protocol))],
        seg("GGG", "GGG", "tso", inferred=inf),
        seg("body", "X" * BODY, placeholder=True),
        seg("B", "B", placeholder=True), seg("polyA", "A" * POLYA),
        seg("GT", revcomp(s2.DT_LINKER)),
        seg("H3", revcomp(H), "tso")],
        name="SMART-seq full-length cDNA before ISPCR")

def _entry(which: str) -> Segment:
    return seg(which, nx.S5 if which == "s5" else nx.S7, which)


def tagmented(a: str, b: str) -> Scene:
    """One Tn5 product, before gap fill.

    Each fragment strand carries a transferred adaptor (entry + ME) on its 5' end. The
    non-transferred strand (ME_RC) stays annealed to each ME, and Tn5's 9-bp stagger leaves
    each fragment 3' end 9 nt short of it: a 9-nt gap opposite the other strand's first 9
    insert bases.
    """
    g = nx.TAGMENTATION_GAP
    me = seg("ME", nx.ME, "me")
    top = [_entry(a), me, seg("ss", "X" * g, placeholder=True),
           seg("core", "XXXXXX...XXXXXX", placeholder=True)]
    bot = [_entry(b), me, seg("ss", "x" * g, placeholder=True),
           *complement_segments([top[-1]])]
    sc = Scene()
    sc.strand("top", top, label="")
    sc.anneal("bottom", bot, to="top", pair=("core'", "core"), label="")
    sc.anneal("nt-top", complement_segments([me], suffix=""), to="top", pair=("ME", "ME"),
              label="")
    sc.anneal("nt-bottom", complement_segments([me], suffix=""), to="bottom",
              pair=("ME", "ME"), label="", above=True)
    return sc


def gap_filled() -> Scene:
    g = nx.TAGMENTATION_GAP
    top = [_entry("s5"), seg("ME", nx.ME, "me"), seg("ss", "X" * g, placeholder=True),
           seg("core", "XXXXXX...XXXXXX", placeholder=True),
           seg("fill", "X" * g, placeholder=True),
           seg("ME_RC", nx.ME_RC, "me"), seg("s7rc", nx.S7_RC, "s7")]
    sc = Scene()
    sc.strand("top", top, label="")
    sc.anneal("bottom", complement_segments(top), to="top", pair=("core'", "core"), label="")
    sc.mark("top", "fill", "top 3' end: gap filled, ME' displaced, copied to the end",
            through="s7rc")
    sc.mark("bottom", "ss'", "bottom 3' end: likewise", through="s5'")
    return sc
