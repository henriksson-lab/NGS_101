#!/usr/bin/env python3
"""Shared drawing functions for the individual SMART-seq protocol pages."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import nextera as nx
import rt
import seqprimers as sp
import smartseq as ss
from chemdraw import (Row, Scene, Segment, annotation_rows, complement_segments, oligo, panel,
                      revcomp)
from illumina import P5, P7
from page import head, info, table

PRE = 4
H = rt.SMART_HANDLE


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


def chunks(*cs):
    return [tuple(list(c) + [None, False])[:3] for c in cs]


# =============================================================================== preamble
PAPERS = {
    "SMART-seq": ("10.1038/nbt.2282", "SMART-seq"),
    "SMART-seq2": ("10.1038/nmeth.2639", "SMART-seq2"),
    "SMART-seq3": ("10.1038/s41587-020-0497-0", "SMART-seq3"),
    "SMART-seq3xpress": ("10.1038/s41587-022-01311-4", "SMART-seq3xpress"),
    "FLASH-seq": ("10.1038/s41587-022-01312-3", "FLASH-seq"),
}


def preamble(protocol: str) -> str:
    doi, title = PAPERS[protocol]
    return f"""<div class="wrap">
<h1>{title} &mdash; full-length single-cell mRNA</h1>

{info(f'''Defining source: <a href="https://doi.org/{doi}">doi:{doi}</a>.''')}

<details class="background">
<summary>Background and interpretation</summary>
<p>Template switching installs the 5' handle when MMLV adds untemplated Cs at the end of
the RNA template. Tn5 then fragments the amplified cDNA and adds the sequencing-adapter
entry points. Cell identity is supplied by the well's i5/i7 index pair rather than by a
cell barcode in the molecule.</p>
</details>
"""


# ====================================================== SMART-seq / SMART-seq2
def ss2_oligos(protocol: str) -> str:
    tail = rt.TSO_G_TAIL_LNA if protocol == "SMART-seq2" else rt.TSO_G_TAIL
    rows = [
        oligo("oligo-dTVN", [seg("", H, "tso"), seg("", ss.SS2_DT_LINKER),
                             seg("", "T" * ss.DT_LEN, None, placeholder=True),
                             seg("", ss.DT_ANCHOR, None, placeholder=True)]),
        oligo("Template Switching Oligo (TSO)",
              [seg("", H, "tso"), seg("", ss.SS2_TSO_LINKER),
               seg("", tail, None, placeholder=True)]),
        oligo("ISPCR primer", [seg("", H, "tso")]),
        oligo("Nextera mosaic end (ME)", [seg("", nx.ME, "me")]),
        oligo("Nextera N/S5xx entry point (s5)", [seg("", nx.S5, "s5")]),
        oligo("Nextera N7xx entry point (s7)", [seg("", nx.S7, "s7")]),
        oligo("Nextera (XT) N/S5xx index primer",
              [seg("", P5, "p5"), seg("", "[8-bp i5]", None, placeholder=True), seg("", nx.S5, "s5")]),
        oligo("Nextera (XT) N7xx index primer",
              [seg("", P7, "p7"), seg("", "[8-bp i7]", None, placeholder=True), seg("", nx.S7, "s7")]),
    ]
    return ("<h2>Adapter and primer sequences</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>\n"
            + info("""The oligo-dT primer, TSO and ISPCR primer share the same 23-nt
            handle. <code>rG</code> denotes RNA.""" +
                   (" <code>+G</code> denotes the SMART-seq2 LNA." if protocol == "SMART-seq2"
                    else "")))


# ------------------------------------------------ step drawings, placed by pairing
# Every duplex panel below is a chemdraw.Scene: strands are given 5'->3' as ordered or
# synthesised and placed by which segment pairs with which. No column is typed by hand,
# and a Scene that would draw an unpaired column raises instead of rendering.

BODY = 20           # drawn transcript body
POLYA = 30          # drawn poly(A): exactly the oligo-dT's (T)30, so the handle overhangs


def _mrna() -> list[Segment]:
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
    s1.strand("mRNA", _mrna())
    s1.anneal("oligo-dT", _dt_primer(handle, link), to="mRNA", pair=("dT", "polyA"))
    s1.mark("oligo-dT", "N", "VN anchor: sits at the poly(A) junction", through="V")
    s1.arrow("oligo-dT", "MMLV")

    s2 = Scene()
    s2.strand("mRNA", _mrna())
    s2.anneal("cDNA", _cdna(handle, link), to="mRNA", pair=("dT", "polyA"))
    s2.mark("cDNA", "CCC", "untemplated CCC, overhanging the mRNA 5' end")

    s3 = Scene()
    s3.strand("mRNA", _mrna())
    s3.anneal("cDNA", _cdna(handle, link), to="mRNA", pair=("dT", "polyA"))
    s3.anneal("TSO", tso, to="cDNA", pair=("rGrG", "CCC"), above=True)
    s3.mark("cDNA", "CCC", "CCC pairs the TSO's G tail")
    s3.arrow("cDNA", "RT switches template")
    return s1, s2, s3


def _ds_cdna(top: list[Segment], cdna: list[Segment], tso_handle: list[Segment],
             fwd: list[Segment], fwd_pair: tuple[str, str], rev_: list[Segment],
             rev_pair: tuple[str, str]) -> list[Row]:
    """Amplified cDNA: bottom = first strand + its copy of the TSO, top = its complement
    written out, both primers placed by the segment they anneal to."""
    sc = Scene()
    sc.strand("top", top, label="")
    sc.anneal("bottom", [*cdna, *complement_segments(tso_handle)], to="top",
              pair=("CCC", "GGG"), label="")
    sc.anneal("fwd", fwd, to="bottom", pair=fwd_pair, label="")
    sc.anneal("rev", rev_, to="top", pair=rev_pair, label="")
    sc.arrow("fwd", "")
    sc.arrow("rev", "")
    sc.stack("fwd", "top", "bottom", "rev")
    return sc.rows()


SS2_HANDLE = seg("H", H, "tso")
SS2_TSO = [SS2_HANDLE, seg("ACAT", ss.SS2_TSO_LINKER), seg("rGrG", "GGG", "tso")]   # rGrG+G, as GGG


def ss2_rt() -> tuple[Scene, ...]:
    return _rt_scenes(SS2_HANDLE, [seg("AC", ss.SS2_DT_LINKER)], SS2_TSO)


def ss2_pcr() -> list[Row]:
    top = [seg("H5", H, "tso"), seg("ACAT", ss.SS2_TSO_LINKER), seg("GGG", "GGG", "tso"),
           seg("body", "X" * BODY, placeholder=True), seg("B", "B", placeholder=True),
           seg("polyA", "A" * POLYA), seg("GT", revcomp(ss.SS2_DT_LINKER)), seg("H3", revcomp(H), "tso")]
    return _ds_cdna(top, _cdna(SS2_HANDLE, [seg("AC", ss.SS2_DT_LINKER)]), SS2_TSO[:2],
                    [seg("ISPCR", H, "tso")], ("ISPCR", "H'"),
                    [seg("ISPCR", H, "tso")], ("ISPCR", "H3"))


def ss2_steps(protocol: str = "SMART-seq2") -> str:
    s1, s2, s3 = ss2_rt()
    p1 = panel(s1.rows(), cls="small",
               caption="(1) The anchored oligo-dTVN anneals at the poly(A) junction; "
                       "MMLV reverse transcribes.")
    p2 = panel(s2.rows(), cls="small",
               caption="(2) Running off the 5' end, MMLV's terminal transferase adds "
                       "three untemplated C.")
    p3 = panel(s3.rows(), cls="small",
               caption=f"(3) The TSO's {'rGrG+G' if protocol == 'SMART-seq2' else 'rGrGrG'} "
                       "pairs with that CCC overhang and the "
                       "polymerase switches template, copying the handle.")
    p4 = panel(ss2_pcr(), cls="long",
               caption="(4) ISPCR amplifies the cDNA from the identical handles at both ends.")
    return "<h2>Step-by-step library generation</h2>\n" + p1 + p2 + p3 + p4 + tagmentation()


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


def tagmentation() -> str:
    rows = []
    for a, b, amp, why in nx.TAGMENTATION_OUTCOMES:
        rows.append(Row(chunks=chunks((f"Product {a}/{b}: {why}",))))
        t = tagmented(a, b)
        t.same_line("top", "nt-bottom")
        t.same_line("bottom", "nt-top")
        rows.extend(t.rows())
        rows.append(Row())
    p5 = panel(rows, cls="long",
               caption="(5) Nextera tagmentation. The Tn5 dimer inserts both adaptors at "
                       "random, leaving a 9-bp gap at each end.")
    p6 = panel(gap_filled().rows(), cls="long",
               caption="(6) 72 C gap fill-in -- the first hold of the Nextera PCR, before any "
                       "denaturation. It is a fill-in, not an extension; skip it and the "
                       "library is lost.")
    return p5 + p6


def ss2_final() -> str:
    con = ss.smartseq2_library()
    duplex = Scene.duplex(con.segments)
    duplex.strands["top"].label = duplex.strands["bottom"].label = ""
    rows = duplex.rows() + annotation_rows(con, prefix_width=PRE)
    return ("<h3>(7&ndash;8) Index PCR, and the final library</h3>\n"
            + panel(rows, cls="long")
            + info(f"""{len(con)}&nbsp;bp excluding the insert; no UMI or cell barcode."""))


# ========================================================== SMART-seq3 family
def ss3_tso_segs(variant: str) -> list[Segment]:
    """The TSO drawn from the same constants as smartseq.ss3_tso()."""
    spacer = ss.SS3_TSO_SPACERS[variant]
    return [seg("", ss.SS3_TSO_ME3, "me"), seg("", ss.SS3_TSO_TAG, "r1"),
            seg("", f"[{ss.SS3_UMI_LEN}-bp UMI]", "umi", placeholder=True),
            *([seg("", spacer, None, placeholder=True)] if spacer else []),
            seg("", rt.TSO_G_TAIL, None, placeholder=True)]


def ss3_oligos(variant: str) -> str:
    tso_names = {"SMART-seq3": "Smartseq3_N8_TSO",
                 "SMART-seq3xpress": "Smartseq3xpress_TSO",
                 "FLASH-seq": "FLASH-seq_TSO"}
    rows = [
        oligo("Smartseq3_OligodT30VN",
              [seg("", ss.SS3_OLIGO_DT_HANDLE, "r3"),
               seg("", "T" * ss.DT_LEN, None, placeholder=True),
               seg("", ss.DT_ANCHOR, None, placeholder=True)],
              mods="/5Biosg/"),
        oligo(tso_names[variant], ss3_tso_segs(variant), mods="/5Biosg/"),
        oligo("Fwd_PCR_primer",
              [seg("", nx.S5, "s5"), seg("", nx.ME, "me"), seg("", ss.SS3_TSO_TAG, "r1")]),
        oligo("Rev_PCR_primer", [seg("", ss.SS3_OLIGO_DT_HANDLE, "r3")]),
    ]
    return ("<h2>Adapter and primer sequences</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>\n"
            + info(f"""The forward and reverse PCR handles differ. The TSO contains the
            terminal 8 nt of the Nextera ME, the 11-bp 5' tag and the 8-bp UMI."""))


SS3_DT_HANDLE = seg("dT handle", ss.SS3_OLIGO_DT_HANDLE, "r3")
SS3_TSO = [seg("ME3", ss.SS3_TSO_ME3, "me"), seg("tag", ss.SS3_TSO_TAG, "r1"),
           seg("UMI", "N" * ss.SS3_UMI_LEN, "umi", placeholder=True),
           seg("rGrG", "GGG", "tso")]                                  # rGrGrG, as GGG


def ss3_tso_scene(variant: str = "SMART-seq3") -> list[Segment]:
    spacer = ss.SS3_TSO_SPACERS[variant]
    return [*SS3_TSO[:3],
            *([seg("spacer", spacer, "r1", placeholder=True)] if spacer else []),
            SS3_TSO[-1]]


def ss3_rt(variant: str = "SMART-seq3") -> tuple[Scene, ...]:
    return _rt_scenes(SS3_DT_HANDLE, [], ss3_tso_scene(variant))


def ss3_pcr(variant: str = "SMART-seq3") -> list[Row]:
    top = [*ss3_tso_scene(variant)[:-1], seg("GGG", "GGG", "tso"),
           seg("body", "X" * BODY, placeholder=True), seg("B", "B", placeholder=True),
           seg("polyA", "A" * POLYA),
           seg("handle3", revcomp(ss.SS3_OLIGO_DT_HANDLE), "r3")]
    fwd = [seg("s5", nx.S5, "s5"), seg("ME", nx.ME, "me"), seg("tag", ss.SS3_TSO_TAG, "r1")]
    return _ds_cdna(top, _cdna(SS3_DT_HANDLE, []), ss3_tso_scene(variant)[:-1],
                    fwd, ("tag", "tag'"),
                    [seg("Rev", ss.SS3_REV_PCR, "r3")], ("Rev", "handle3"))


def ss3_steps(variant: str = "SMART-seq3") -> str:
    s1, s2, s3 = ss3_rt(variant)
    p1 = panel(s1.rows(), cls="small",
               caption="(1) The anchored oligo-dT primer anneals at the poly(A) junction; "
                       "MMLV reverse transcribes.")
    p2 = panel(s2.rows(), cls="small",
               caption="(2) At the mRNA 5' end, MMLV adds three untemplated C.")
    p3 = panel(s3.rows(), cls="small",
               caption="(3) Template switching, but the TSO now carries an 11-bp tag and an "
                       "8-bp UMI. Both are attached BEFORE any amplification, and only at a "
                       "genuine transcript 5' end.")
    p4 = panel(ss3_pcr(variant), cls="long",
               caption="(4) Two different primers. Amplification is no longer single-primer, "
                       "and no longer suppressive. The forward primer's s5 + ME head overhangs: "
                       "only its last 8 nt of ME and the tag anneal.")
    return ("<h2>Step-by-step library generation</h2>\n"
            + p1 + p2 + p3 + p4 + tagmentation())


def ss3_final(variant: str = "SMART-seq3") -> str:
    five = ss.smartseq3_library(variant)
    internal = ss.smartseq3_library(five_prime=False)
    five_duplex = Scene.duplex(five.segments)
    five_duplex.strands["top"].label = five_duplex.strands["bottom"].label = ""
    internal_duplex = Scene.duplex(internal.segments)
    internal_duplex.strands["top"].label = internal_duplex.strands["bottom"].label = ""
    r5 = five_duplex.rows() + annotation_rows(five, prefix_width=PRE)
    ri = internal_duplex.rows()
    return ("<h3>(9) Final library structures</h3>\n"
            + "<h3>5' fragments &mdash; these retained the TSO</h3>\n" + panel(r5, cls="long")
            + info("""Contains the <r1>11-bp 5' tag</r1> followed by the
            <umi>8-bp UMI</umi>.""")
            + "<h3>Internal fragments &mdash; these did not</h3>\n" + panel(ri, cls="long"))


# =============================================================== sequencing
def seq_primer(primer: list[Segment], on: str, pair: tuple[str, str]) -> list[Row]:
    """The SMART-seq2 library duplex with one sequencing primer annealed to strand `on`."""
    con = ss.smartseq2_library()
    sc = Scene()
    sc.strand("top", list(con), label="")
    sc.anneal("bottom", complement_segments(list(con)), to="top", pair=("cDNA'", "cDNA"),
              label="")
    sc.anneal("primer", primer, to=on, pair=pair, label="")
    sc.arrow("primer", "")
    sc.stack(*(("primer", "top", "bottom") if on == "bottom" else ("top", "bottom", "primer")))
    return sc.rows()


def seq_primer_drawings() -> dict[str, tuple[list[Segment], str, tuple[str, str]]]:
    """Segment drawings of the four stock Nextera primers: (segments, strand, pair).
    Built from the same nextera.py pieces; the self-test asserts each spells exactly the
    sp.NEXTERA primer of the same key."""
    s5, me, s7 = seg("s5", nx.S5, "s5"), seg("ME", nx.ME, "me"), seg("s7", nx.S7, "s7")
    me_rc = seg("ME_RC", nx.ME_RC, "me")
    return {"R1": ([s5, me], "bottom", ("s5", "s5'")),
            "I1": ([me_rc, seg("s7rc", nx.S7_RC, "s7")], "bottom", ("ME_RC", "ME '")),
            "I2": ([me_rc, seg("s5rc", nx.S5_RC, "s5")], "top", ("ME_RC", "ME")),
            "R2": ([s7, me], "top", ("s7", "s7"))}


def sequencing(protocol: str) -> str:
    if protocol in ("SMART-seq", "SMART-seq2"):
        tables = sp.section(ss.smartseq2_library(), ss.SS2_SEQ_PRIMERS,
                            heading=f"Sequencing primers: {protocol}")
    else:
        tables = (sp.section(ss.smartseq3_library(protocol), ss.SS3_SEQ_PRIMERS,
                             heading=f"Sequencing primers: {protocol} 5' fragments")
                  + sp.section(ss.smartseq3_library(protocol, five_prime=False),
                               ss.SS3_SEQ_PRIMERS,
                               heading=f"Sequencing primers: {protocol} internal fragments"))
    return ('<h2 id="seq-primers">Library sequencing</h2>\n'
            + info("""All four reads use the stock Nextera sequencing primers. Primer
            positions and the first bases read are shown below.""")
            + tables
            + info("""The i5 + i7 index pair identifies the source well.""")
            + sources(protocol)
            + "</div>")


def sources(protocol: str) -> str:
    related = {
        "SMART-seq": [],
        "SMART-seq2": [("SMART-seq", PAPERS["SMART-seq"][0])],
        "SMART-seq3": [("SMART-seq2", PAPERS["SMART-seq2"][0])],
        "SMART-seq3xpress": [("SMART-seq3", PAPERS["SMART-seq3"][0])],
        "FLASH-seq": [("SMART-seq3", PAPERS["SMART-seq3"][0])],
    }[protocol]
    links = ''.join(f'<li><a href="https://doi.org/{doi}">{name}</a></li>'
                    for name, doi in related)
    links += ('<li><a href="https://teichlab.github.io/scg_lib_structs/methods_html/'
              'SMART-seq_family.html">scg_lib_structs schematic</a></li>')
    return ('<details class="sources"><summary>Related sources</summary><ul>'
            + links + '</ul></details>')


def render_page(protocol: str) -> str:
    if protocol in ("SMART-seq", "SMART-seq2"):
        body = [ss2_oligos(protocol), ss2_steps(protocol), ss2_final()]
    else:
        body = [ss3_oligos(protocol), ss3_steps(protocol), ss3_final(protocol)]
    return "\n".join([head(f"{protocol} library chemistry"), preamble(protocol),
                      *body, sequencing(protocol)])
