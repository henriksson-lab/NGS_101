#!/usr/bin/env python3
"""Generate smartseq.html -- the SMART-seq family chemistry page."""
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
                      revcomp, strand_row)
from illumina import P5, P7
from page import caveat, head, info, legend, table

OUT = HERE.parent / "smartseq.html"
PRE = 4
H = rt.SMART_HANDLE


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


def chunks(*cs):
    return [tuple(list(c) + [None, False])[:3] for c in cs]


# =============================================================================== preamble
def preamble() -> str:
    return f"""<div class="wrap">
<h1>SMART-seq family &mdash; full-length single-cell mRNA</h1>

{info("""Five related protocols that all read the <b>whole transcript</b> rather than
counting its 3' end: <b>SMART-seq</b>, <b>SMART-seq2</b>, <b>SMART-seq3</b>,
<b>SMART-seq3xpress</b> and <b>FLASH-seq</b>. Structures transcribed from
<a href="https://teichlab.github.io/scg_lib_structs/methods_html/SMART-seq_family.html">
scg_lib_structs</a>.""")}

{info("""They divide into two architectures. <b>SMART-seq / SMART-seq2</b> put the
<i>same</i> handle on both ends of the cDNA, so one primer amplifies it. <b>SMART-seq3</b>
and its descendants deliberately make the two ends <i>different</i> and hide an 11-bp tag
and a UMI in the template-switching oligo &mdash; which is what lets a single protocol give
full-length coverage <i>and</i> molecule counting.""")}

{legend("""<b>Two ideas carry this whole page.</b>
<b>Template switching</b> puts a defined handle on the 5' end of a cDNA with no ligation,
by exploiting the untemplated Cs that MMLV adds when it runs off the end of the template.
<b>Tn5 tagmentation</b> then fragments the amplified cDNA and installs the sequencing
adapters in one step. Both are reused across many protocols and are written up separately
in <code>ref/concepts/</code>.""")}

{caveat("""<b>Neither chemistry carries a cell barcode.</b> One cell goes in one well, and
the cell's identity is the <i>i5 + i7 index pair</i> of its library. That is what caps these
methods at plate scale &mdash; and what the split-pool methods exist to fix.""")}
"""


# ====================================================== SMART-seq / SMART-seq2
def ss2_oligos() -> str:
    rows = [
        oligo("oligo-dTVN", [seg("", H, "tso"), seg("", ss.SS2_DT_LINKER),
                             seg("", "T" * ss.DT_LEN, None, placeholder=True),
                             seg("", ss.DT_ANCHOR, None, placeholder=True)]),
        oligo("Template Switching Oligo (TSO)",
              [seg("", H, "tso"), seg("", ss.SS2_TSO_LINKER),
               seg("", rt.TSO_G_TAIL_LNA, None, placeholder=True)]),
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
            + info("""<b>The same 23-nt handle appears three times</b> &mdash; on the
            oligo-dT primer, on the TSO, and as the PCR primer. That is the entire design:
            after template switching both ends of the cDNA are identical, so a single
            primer amplifies it. <code>+G</code> is a locked nucleic acid, the SMART-seq2
            improvement; <code>rG</code> are ribonucleotides. The four sequencing primers
            are listed under <a href="#seq-primers">Library sequencing</a>."""))


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


def ss2_steps() -> str:
    s1, s2, s3 = ss2_rt()
    p1 = panel(s1.rows(), cls="small",
               caption="(1) The anchored oligo-dTVN anneals at the poly(A) junction; "
                       "MMLV reverse transcribes.")
    p2 = panel(s2.rows(), cls="small",
               caption="(2) Running off the 5' end, MMLV's terminal transferase adds "
                       "three untemplated C.")
    p3 = panel(s3.rows(), cls="small",
               caption="(3) The TSO's rGrG+G pairs with that CCC overhang and the "
                       "polymerase switches template, copying the handle.")
    p4 = panel(ss2_pcr(), cls="long",
               caption="(4) One primer, both ends. Identical ends make this "
                       "semi-suppressive PCR: short products form a pan-handle that "
                       "competes with priming, biasing amplification toward full-length cDNA.")
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
    return p5 + caveat("""<b>Only one product in three survives.</b> A fragment with
    <s5>s5</s5> at both ends has P5 at both ends and no P7 site; <s7>s7</s7> at both ends has
    the converse. Only the <s5>s5</s5>/<s7>s7</s7> heteroduplex can bridge P5 to P7, so the
    other two amplify linearly at best. This is suppression by design.""") + p6


def ss2_final() -> str:
    con = ss.smartseq2_library()
    rows = [strand_row(con, "top"), strand_row(con, "bottom")] + annotation_rows(con, prefix_width=PRE)
    return ("<h3>(7&ndash;8) Index PCR, and the final library</h3>\n"
            + panel(rows, cls="long")
            + info(f"""{len(con)}&nbsp;bp excluding the insert. <b>No UMI, no cell
            barcode.</b> Reads land anywhere in the transcript, in either orientation &mdash;
            there is no positional information anywhere in this construct."""))


# ========================================================== SMART-seq3 family
def ss3_tso_segs(variant: str) -> list[Segment]:
    """The TSO drawn from the same constants as smartseq.ss3_tso()."""
    spacer = ss.SS3_TSO_SPACERS[variant]
    return [seg("", ss.SS3_TSO_ME3, "me"), seg("", ss.SS3_TSO_TAG, "r1"),
            seg("", f"[{ss.SS3_UMI_LEN}-bp UMI]", "umi", placeholder=True),
            *([seg("", spacer, None, placeholder=True)] if spacer else []),
            seg("", rt.TSO_G_TAIL, None, placeholder=True)]


def ss3_oligos() -> str:
    rows = [
        oligo("Smartseq3_OligodT30VN",
              [seg("", ss.SS3_OLIGO_DT_HANDLE, "r3"),
               seg("", "T" * ss.DT_LEN, None, placeholder=True),
               seg("", ss.DT_ANCHOR, None, placeholder=True)],
              mods="/5Biosg/"),
        *[oligo(name, ss3_tso_segs(v), mods="/5Biosg/")
          for name, v in (("Smartseq3_N8_TSO", "SMART-seq3"),
                          ("Smartseq3xpress_TSO", "SMART-seq3xpress"),
                          ("FLASH-seq_TSO", "FLASH-seq"))],
        oligo("Fwd_PCR_primer",
              [seg("", nx.S5, "s5"), seg("", nx.ME, "me"), seg("", ss.SS3_TSO_TAG, "r1")]),
        oligo("Rev_PCR_primer", [seg("", ss.SS3_OLIGO_DT_HANDLE, "r3")]),
    ]
    return ("<h2>Adapter and primer sequences</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>\n"
            + info(f"""The two handles now <b>differ</b>, so PCR uses two primers rather than
            one. Note what the TSO actually is: <me>{ss.SS3_TSO_ME3}</me> is the 3' end of the
            Nextera <me>mosaic end</me>, so the forward PCR primer
            (<s5>s5</s5>&nbsp;+&nbsp;<me>ME</me>&nbsp;+&nbsp;<r1>tag</r1>) rebuilds a
            complete Nextera s5 arm as it amplifies &mdash; the TSO was already most of one.
            The three family members differ <i>only</i> in the spacer between UMI and G tail:
            {", ".join(f"<code>{x}</code>" if x else "none"
                       for x in ss.SS3_TSO_SPACERS.values())}."""))


SS3_DT_HANDLE = seg("dT handle", ss.SS3_OLIGO_DT_HANDLE, "r3")
SS3_TSO = [seg("ME3", ss.SS3_TSO_ME3, "me"), seg("tag", ss.SS3_TSO_TAG, "r1"),
           seg("UMI", "N" * ss.SS3_UMI_LEN, "umi", placeholder=True),
           seg("rGrG", "GGG", "tso")]                                  # rGrGrG, as GGG


def ss3_rt() -> tuple[Scene, ...]:
    return _rt_scenes(SS3_DT_HANDLE, [], SS3_TSO)


def ss3_pcr() -> list[Row]:
    top = [*SS3_TSO[:3], seg("GGG", "GGG", "tso"),
           seg("body", "X" * BODY, placeholder=True), seg("B", "B", placeholder=True),
           seg("polyA", "A" * POLYA),
           seg("handle3", revcomp(ss.SS3_OLIGO_DT_HANDLE), "r3")]
    fwd = [seg("s5", nx.S5, "s5"), seg("ME", nx.ME, "me"), seg("tag", ss.SS3_TSO_TAG, "r1")]
    return _ds_cdna(top, _cdna(SS3_DT_HANDLE, []), SS3_TSO[:3],
                    fwd, ("tag", "tag'"),
                    [seg("Rev", ss.SS3_REV_PCR, "r3")], ("Rev", "handle3"))


def ss3_steps() -> str:
    _, _, s3 = ss3_rt()
    p3 = panel(s3.rows(), cls="small",
               caption="(3) Template switching, but the TSO now carries an 11-bp tag and an "
                       "8-bp UMI. Both are attached BEFORE any amplification, and only at a "
                       "genuine transcript 5' end.")
    p4 = panel(ss3_pcr(), cls="long",
               caption="(4) Two different primers. Amplification is no longer single-primer, "
                       "and no longer suppressive. The forward primer's s5 + ME head overhangs: "
                       "only its last 8 nt of ME and the tag anneal.")
    return "<h2>Step-by-step library generation</h2>\n" + p3 + p4 + info(
        """Steps 1, 2 and 5&ndash;8 are as for SMART-seq2: oligo-dT priming and reverse
        transcription (here with <b>Maxima H&minus;</b>, in NaCl with 5% PEG), the untemplated
        CCC, then purification, Nextera tagmentation, the 72&nbsp;&deg;C gap fill-in and the
        index PCR.""")


def ss3_final() -> str:
    five = ss.smartseq3_library()
    internal = ss.smartseq3_library(five_prime=False)
    r5 = [strand_row(five, "top"), strand_row(five, "bottom")] + annotation_rows(five, prefix_width=PRE)
    ri = [strand_row(internal, "top"), strand_row(internal, "bottom")]
    variants = table(
        ["Method", "TSO spacer", "5' fragment length vs SMART-seq2"],
        [[v, f"<code>{s or '&mdash;'}</code>",
          f"+{len(ss.SS3_TSO_TAG) + ss.SS3_UMI_LEN + 3 + len(s)}&nbsp;bp"]
         for v, s in ss.SS3_TSO_SPACERS.items()])
    return ("<h3>(9) Final library structures</h3>\n"
            + "<h3>5' fragments &mdash; these retained the TSO</h3>\n" + panel(r5, cls="long")
            + info("""Carries the <r1>11-bp tag</r1> and the <umi>UMI</umi>. A read beginning
            with the tag is known to come from the transcript's 5' end, and the following
            8&nbsp;nt count the molecule it came from.""")
            + "<h3>Internal fragments &mdash; these did not</h3>\n" + panel(ri, cls="long")
            + caveat("""<b>This is character-for-character a SMART-seq2 fragment</b> &mdash;
            no tag, no UMI, no positional information. The page's self-test asserts it.
            Which is exactly why the 11-bp tag has to exist: without it there is no way to
            tell a 5' fragment from an internal one, and therefore no way to count
            molecules.""")
            + "<h3>The three variants</h3>\n" + variants)


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


def sequencing() -> str:
    d = seq_primer_drawings()
    for k, (segs, _, _) in d.items():
        if "".join(x.top for x in segs) != sp.NEXTERA[k].seq:
            raise ValueError(f"{k} drawing does not spell sp.NEXTERA[{k!r}]")
    p = {k: seq_primer(*v) for k, v in d.items()}
    tables = (sp.section(ss.smartseq2_library(), ss.SS2_SEQ_PRIMERS,
                         heading="Sequencing primers: SMART-seq2 (and SMART-seq3 internal fragments)")
              + sp.section(ss.smartseq3_library(), ss.SS3_SEQ_PRIMERS,
                           heading="Sequencing primers: SMART-seq3 5' fragments"))
    return ('<h2 id="seq-primers">Library sequencing</h2>\n'
            + info("""Stock Nextera: all four primers are built from the same two pieces,
            the <me>mosaic end</me> and an entry point. The tables below are computed from
            the final libraries; a SMART-seq3 <i>internal</i> fragment is identical to a
            SMART-seq2 fragment, and SMART-seq3xpress / FLASH-seq 5' fragments read the
            same first bases as SMART-seq3 (their spacer follows the UMI). The self-test
            verifies the primers on every one of them. Note Read&nbsp;1 on a SMART-seq3
            5' fragment: the TSO supplied the last 8&nbsp;nt of the <me>ME</me>, so the
            stock primer ends exactly where the <r1>tag</r1> begins.""")
            + tables
            + "<h3>(1) Read 1 &mdash; bottom strand as template</h3>\n"
            + panel(p["R1"], cls="long")
            + "<h3>(2) Index 1 (i7) &mdash; bottom strand as template</h3>\n"
            + panel(p["I1"], cls="long")
            + "<h3>(3) Index 2 (i5) &mdash; top strand as template, after cluster regeneration</h3>\n"
            + panel(p["I2"], cls="long")
            + "<h3>(4) Read 2 &mdash; top strand as template</h3>\n"
            + panel(p["R2"], cls="long")
            + info("""The cell is identified by the <b>i5 + i7 pair</b>, read in steps 2 and
            3 &mdash; one well, one index combination, one cell.""")
            + "</div>")


def main() -> None:
    parts = [head("SMART-seq Family Chemistry"), preamble(),
             '<h1 id="smart-seq">SMART-seq / SMART-seq2</h1>',
             ss2_oligos(), ss2_steps(), ss2_final(),
             '<h1 id="smart-seq3">SMART-seq3 / SMART-seq3xpress / FLASH-seq</h1>',
             ss3_oligos(), ss3_steps(), ss3_final(),
             sequencing()]
    OUT.write_text("\n".join(parts), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
