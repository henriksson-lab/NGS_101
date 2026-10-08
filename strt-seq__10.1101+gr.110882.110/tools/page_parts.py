#!/usr/bin/env python3
"""Shared drawing implementation for the three individual STRT pages."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import strt as S
from chemdraw import (Scene, Segment, annotation_rows, complement_segments, oligo, panel,
                      revcomp, strand_row)
import illumina as il
import nextera as nx
from page import head, info, table
import seqprimers as sp


PAPERS = {
    S.STRT: ("10.1101/gr.110882.110", "STRT-seq"),
    S.C1: ("10.1038/nmeth.2772", "STRT-C1"),
    S.TWO_I: ("10.1038/s41598-017-16546-4", "STRT-seq-2i"),
}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def preamble(protocol: str) -> str:
    doi, title = PAPERS[protocol]
    notes = {
        S.STRT: "01_strt-seq.html",
        S.C1: "../strt-seq__10.1101+gr.110882.110/01_strt-seq.html",
        S.TWO_I: "../strt-seq__10.1101+gr.110882.110/01_strt-seq.html",
    }[protocol]
    caveat = ""
    if protocol == S.STRT:
        caveat = '<div class="caveat"><b>Primer sequence unavailable.</b> The defining paper states that Read 1 uses a custom primer but does not print it. Its landing site is reconstructed from the published final construct.</div>'
    elif protocol == S.C1:
        caveat = '<div class="caveat"><b>Partial methods.</b> The accessible Supplementary Note defines the molecule and says it is sequenced without library PCR; the paywalled Online Methods contain the unavailable selection details.</div>'
    elif protocol == S.TWO_I:
        caveat = '<div class="caveat"><b>Inferred junction.</b> The published P2 PCR primer overlaps the first three bases of the subarray barcode. The dotted barcode/P7 end shows the full-barcode outcome; the paper does not state whether proofreading instead trims the primer or overwrites those bases.</div>'
    return f'''<div class="wrap">
<h1>{title} &mdash; 5'-end single-cell RNA sequencing</h1>
<p class="research-notes"><a href="{notes}">Research notes</a></p>
{info(f'Defining source: <a href="https://doi.org/{doi}">doi:{doi}</a>.')}
{caveat}
'''


def oligos(protocol: str) -> str:
    if protocol == S.STRT:
        rows = [
            oligo("STRT-V3-T30", [seg("SMART handle", S.H, "tso"), seg("SalI junction", S.STRT_DT_LINK), seg("T30", "T" * S.STRT_DT_LEN), seg("VN", "VN", placeholder=True)], mods="biotin"),
            oligo("STRT-V2-n TSO", [seg("SMART handle", S.H, "tso"), seg("spacer", S.STRT_TSO_LINK, "r1"), seg("cell barcode", "B" * 6, "cbc", placeholder=True), seg("rGrGrG", "rGrGrG", placeholder=True)]),
            oligo("STRT-PCR", [seg("SMART handle", S.H, "tso")], mods="biotin"),
            oligo("P2 adapter, top", [seg("P7 fragment", il.P7[:22], "p7"), seg("stem + T", il.TRUSEQ_READ2[-12:], "r2")]),
            oligo("Library PCR 1", [seg("P5 fragment", il.P5[:25], "p5"), seg("SMART handle", S.H, "tso")]),
            oligo("Library PCR 2", [seg("P7 fragment", il.P7[:22], "p7")]),
        ]
    elif protocol == S.C1:
        rows = [
            oligo("C1-P1-T31", [seg("C1 handle", S.C1_HANDLE, "p5"), seg("PvuI junction", S.C1_DT_LINK), seg("T31", "T" * S.C1_DT_LEN)], mods="Bio"),
            oligo("C1-P1-RNA-TSO", [seg("C1 handle", S.C1_HANDLE, "p5"), seg("UMI", "N" * 5, "umi", placeholder=True), seg("rGrGrG", "rGrGrG", placeholder=True)], mods="Bio; all-ribo"),
            oligo("C1-P1-PCR-2", [seg("extra G", "G"), seg("C1 handle", S.C1_HANDLE, "p5")], mods="Bio"),
            oligo("C1-TN5-x", S.tn5_indexed_top()),
            oligo("C1-TN5-U", [seg("paired tail'", S.TN5_U, "me")], mods="5Phos"),
        ]
    else:
        rows = [
            oligo("STRT-P1-T31", [seg("C1 handle", S.C1_HANDLE, "p5"), seg("CG", S.C1_DT_LINK), seg("T31", "T" * S.C1_DT_LEN)], mods="Bio"),
            oligo("P1B-UMI-RNA-TSO", [seg("P1B", S.P1B, "r1"), seg("UMI", "N" * 6, "umi", placeholder=True), seg("rGrGrG", "rGrGrG", placeholder=True)], mods="Bio; all-ribo"),
            oligo("DI-P1A-idx-P1B", [seg("P5", il.P5, "p5"), seg("well index", "W" * 5, "cbc", placeholder=True), seg("P1B minus terminal T", S.P1B[:-1], "r1")], mods="Bio"),
            oligo("STRT-Tn5-Idx", S.tn5_indexed_top()),
            oligo("STRT-TN5-U", [seg("paired tail'", S.TN5_U, "me")], mods="5Phos"),
            oligo("P1/P2 second-PCR primers", [seg("P5 fragment", il.P5[:24], "p5"), seg(" / ", "", placeholder=True), seg("P7", il.P7, "p7")]),
        ]
    return "<h2>Key oligos</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>"


def rt_scenes(protocol: str) -> tuple[Scene, Scene]:
    body = seg("body", "XXXXXXXX...XXXXXXXX", placeholder=True)
    if protocol == S.STRT:
        poly = "A" * S.STRT_DT_LEN
        primer = [seg("handle", S.H, "tso"), seg("CGAC", S.STRT_DT_LINK),
                  seg("dT", "T" * S.STRT_DT_LEN)]
        tso = [seg("TSO handle", S.H, "tso"), seg("spacer", S.STRT_TSO_LINK, "r1"), seg("cell barcode", "B" * 6, "cbc", placeholder=True), seg("GGG", "GGG", "tso")]
    elif protocol == S.C1:
        poly = "A" * S.C1_DT_LEN
        primer = [seg("handle", S.C1_HANDLE, "p5"), seg("CG", S.C1_DT_LINK), seg("dT", "T" * S.C1_DT_LEN)]
        tso = [seg("TSO handle", S.C1_HANDLE, "p5"), seg("UMI", "N" * 5, "umi", placeholder=True), seg("GGG", "GGG", "tso")]
    else:
        poly = "A" * S.C1_DT_LEN
        primer = [seg("handle", S.C1_HANDLE, "p5"), seg("CG", S.C1_DT_LINK), seg("dT", "T" * S.C1_DT_LEN)]
        tso = [seg("P1B", S.P1B, "r1"), seg("UMI", "N" * 6, "umi", placeholder=True), seg("GGG", "GGG", "tso")]
    mrna = [body, seg("poly(A)", poly)]
    first = [*primer, *complement_segments([body]), seg("CCC", "CCC", "tso")]
    p = Scene(); p.strand("mRNA", mrna); p.anneal("oligo-dT", primer, to="mRNA", pair=("dT", "poly(A)"), mod5="biotin"); p.arrow("oligo-dT", "reverse transcriptase")
    sw = Scene(); sw.strand("mRNA", mrna); sw.anneal("first-strand cDNA", first, to="mRNA", pair=("dT", "poly(A)"), mod5="biotin"); sw.anneal("TSO", tso, to="first-strand cDNA", pair=("GGG", "CCC"), above=True); sw.arrow("first-strand cDNA", "template switch")
    return p, sw


def indexed_tn5_scene() -> Scene:
    top = S.tn5_indexed_top()
    sc = Scene(); sc.strand("indexed top", top, label="indexed strand"); sc.anneal("TN5-U", [seg("paired tail'", S.TN5_U, "me")], to="indexed top", pair=("paired tail'", "GCGTC"), label="TN5-U", mod5="p")
    sc.mark("indexed top", "subarray barcode", "8-nt Tn5 barcode")
    return sc


def strt_p2_scene() -> Scene:
    """The original STRT TSO-end fragment after P2-adapter ligation."""
    top = [seg("TSO end", S.H, "tso"),
           seg("cell barcode", "B" * 6, "cbc", placeholder=True),
           seg("5′ cDNA fragment", "XXXXXXXX...XXXXXXXX", placeholder=True),
           seg("P2 adapter'", revcomp(S.P2_TOP), "p7")]
    sc = Scene.duplex(top)
    sc.junction("top", "5′ cDNA fragment", "P2 adapter'", "P2 ligation")
    return sc


def final_panel(protocol: str) -> str:
    lib = S.library(protocol)
    caption = "Final Read-1-sense library molecule."
    if protocol == S.TWO_I:
        caption = "INFERRED — final Read-1-sense molecule under the full-subarray-barcode outcome; dotted regions mark the unresolved P2-primer junction."
    return panel([strand_row(lib, "top"), strand_row(lib, "bottom"), *annotation_rows(lib)], caption=caption)


def steps(protocol: str) -> str:
    p, sw = rt_scenes(protocol)
    common = f'''<h2>Library generation</h2>
<h3>(1) Reverse transcription</h3>{panel(p.rows(), cls="small", caption="The biotinylated oligo-dT primes the poly(A) tail.")}
<h3>(2) Template switching at the RNA 5' end</h3>{panel(sw.rows(), cls="small", caption="The TSO's terminal GGG pairs with the cDNA's untemplated CCC, placing the cell barcode or UMI immediately before the transcript 5' end.")}'''
    if protocol == S.STRT:
        mid = '''<h3>(3) Pool wells and amplify from the shared SMART handles</h3>
''' + panel(Scene.duplex([seg("SMART handle / TSO end", S.H, "tso"), seg("cell barcode", "B" * 6, "cbc", placeholder=True), seg("5' cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True), seg("SMART handle / oligo-dT end", S.H, "tso")]).rows(), caption="Biotinylated single-primer PCR preserves the TSO barcode; all 96 wells can now be handled together.") + '''
<h3>(4) Immobilise, fragment, repair, A-tail and ligate the P2 adapter</h3>
''' + panel(strt_p2_scene().rows(), caption="SalI releases oligo-dT-end fragments at the published handle/poly(T) site; bead-bound TSO-end fragments receive P2 and remain for library PCR. ** marks the adapter ligation.")
    elif protocol == S.C1:
        mid = '''<h3>(3) Amplify full-length cDNA from the shared C1 handles</h3>
''' + panel(Scene.duplex([seg("C1 handle / TSO end", S.C1_HANDLE, "p5"), seg("UMI", "N" * 5, "umi", placeholder=True), seg("5' cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True), seg("C1 handle / oligo-dT end", S.C1_HANDLE, "p5")]).rows(), caption="The all-RNA TSO is not a competing PCR primer; C1-P1-PCR-2 supplies biotin.") + '''
<h3>(4) Tagment with one barcoded adapter per C1 chamber</h3>
''' + panel(indexed_tn5_scene().rows(), cls="small", caption="The adapter contributes P7, the 8-nt cell barcode and the Read-2/ME end.") + '''
<h3>(5) Select the TSO-end strand and sequence directly</h3>
''' + panel([strand_row(S.c1_library(), "top"), *annotation_rows(S.c1_library())], caption="No library PCR: the selected single strand already carries the C1 Read-1 end and barcoded Tn5 end.")
    else:
        mid = '''<h3>(3) Add the 5-nt well index during cDNA PCR</h3>
''' + panel(Scene.duplex([seg("P5", il.P5, "p5"), seg("well index", "W" * 5, "cbc", placeholder=True), seg("P1B", S.P1B, "r1"), seg("UMI", "N" * 6, "umi", placeholder=True), seg("5' cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True), seg("P1A end", S.C1_HANDLE, "p5")]).rows(), caption="The indexed primer acts at the TSO end; after its tail is copied, the shared P1A primer amplifies the cDNA.") + '''
<h3>(4) Pool by subarray and tagment with an 8-nt subarray index</h3>
''' + panel(indexed_tn5_scene().rows(), cls="small", caption="Each of 96 subarray pools receives a different indexed Tn5 adapter.") + '''
<h3>(5) Streptavidin capture, single-strand elution and eight-cycle PCR</h3>
''' + panel([strand_row(S.two_i_library(), "top"), *annotation_rows(S.two_i_library())], caption="INFERRED — selected TSO-end strand after the final P1/P2 PCR; dotted regions mark the unresolved P2/barcode junction.")
    final_step = 5 if protocol == S.STRT else 6
    return common + mid + f'<h3>({final_step}) Final library</h3>' + final_panel(protocol)


def sequencing(protocol: str) -> str:
    lib = S.library(protocol)
    if protocol == S.STRT:
        primers = (S.STRT_READ1,)
        hit = sp.locate(lib, S.STRT_READ1)
        rows = [("Read 1", "~55", "6-nt cell barcode &middot; GGG &middot; transcript 5' end")]
        note = "The paper does not publish the custom primer sequence; the drawn landing site ends immediately before the barcode."
    elif protocol == S.C1:
        primers = S.C1_PRIMERS
        hits = [sp.locate(lib, p) for p in S.C1_PRIMERS]
        if any(h is None for h in hits): raise ValueError("C1 read site missing")
        rows = [("Read 1", "50", "5-nt UMI &middot; GGG &middot; transcript 5' end"), ("Index 1", "8", "Tn5 cell barcode")]
        note = "Read 1 and index-primer sequences are not printed in the accessible paper materials; their landing sites follow from the published molecule."
    else:
        primers = S.TWO_I_INDEX_PRIMERS
        hits = [sp.locate(lib, p) for p in S.TWO_I_INDEX_PRIMERS]
        if any(h is None for h in hits): raise ValueError("2i index site missing")
        rows = [("Read 1", "45 or 48", "6-nt UMI &middot; GGG &middot; transcript 5' end"), ("Index 1", "8", "subarray index"), ("Index 2", "5", "well index")]
        note = "The paper reports 45 cycles in Methods and 48 in its supplement. Its printed Read-1 primer contains six N opposite a five-base well index; the constant 3' P1B portion still fixes the read start."
        # Scene enforces the printed degenerate primer against the published P5+5-index+P1B template.
        sc = Scene.duplex([seg("P5", il.P5, "p5"),
                           seg("well index", "W" * 5, "cbc", placeholder=True),
                           seg("P1B", S.P1B, "r1"),
                           seg("UMI", "N" * 6, "umi", placeholder=True)])
        sc.anneal("DI-Read1-Seq", [seg("P1B primer (paired 3' end)", S.P1B, "r1")],
                  to="bottom", pair=("P1B primer (paired 3' end)", "P1B'"))
        sc.arrow("DI-Read1-Seq", "UMI, then transcript")
        note += panel(sc.rows(), cls="small", caption="The printed 6-N Read-1 primer spans the last P5 base plus the five-base well index, then pairs exactly over P1B.")
    if protocol == S.STRT and hit is None:
        raise ValueError("STRT Read 1 site missing")
    return ('<h2>Read layout</h2>' + sp.diagram(lib, primers)
            + table(("Read", "Cycles", "Content"), rows) + info(note) + '</div>')


def render_page(protocol: str, out: Path) -> None:
    title = PAPERS[protocol][1]
    out.write_text("\n".join([head(f"{title} library chemistry"), preamble(protocol), oligos(protocol), steps(protocol), sequencing(protocol)]), encoding="utf-8")
    print(f"wrote {out}  ({out.stat().st_size:,} bytes)")
