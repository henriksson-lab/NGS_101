#!/usr/bin/env python3
"""Shared renderer for individual sci-ATAC-seq protocols."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import sciatac as S
from chemdraw import (Scene, Segment, annotation_rows, complement_segments, oligo, panel,
                      revcomp, strand_row)
import illumina as il
import nextera as nx
from page import head, info, table
import seqprimers as sp

PAPERS = {
    S.SCI15: ("10.1126/science.aab1601", "sci-ATAC-seq"),
    S.SCI18: ("10.1038/nature25981", "sci-ATAC-seq (2018)"),
    S.SCI3: ("10.1126/science.aba7612", "sci-ATAC-seq3"),
}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def preamble(protocol: str) -> str:
    doi, title = PAPERS[protocol]
    if protocol == S.SCI15:
        caveat = '<div class="caveat"><b>Inferred PCR indices.</b> The defining paper’s gated supplement contains the PCR-primer list. Dotted 8-nt PCR indices use the compatible Amini 2014 Nextera-kit layout; their exact 2015 length and bases remain unavailable.</div>'
    elif protocol == S.SCI3:
        caveat = ""
    else:
        caveat = ""
    return f'''<div class="wrap"><h1>{title} &mdash; combinatorial single-cell chromatin accessibility</h1>
{info(f'Defining source: <a href="https://doi.org/{doi}">doi:{doi}</a>.')}{caveat}'''


def two_oligos(protocol: str) -> str:
    pcr_n = 8 if protocol == S.SCI15 else 10
    inf = protocol == S.SCI15
    rows = [
        oligo("T5 indexed transposon", S.t5()), oligo("T7 indexed transposon", S.t7()),
        oligo("pMENTS", [seg("ME'", nx.ME_RC, "me")], mods="5Phos"),
        oligo("Read 1", [seg("extension A", S.R1_EXT, "r1"), seg("ME", nx.ME, "me")]),
        oligo("Read 2", [seg("extension B", S.R2_EXT, "r1"), seg("ME", nx.ME, "me")]),
        oligo("P5 PCR primer", [seg("P5", il.P5, "p5"), seg("i5", "N" * pcr_n, "cbc", placeholder=True, inferred=inf), seg("s5", nx.S5, "s5")]),
        oligo("P7 PCR primer", [seg("P7", il.P7, "p7"), seg("i7", "N" * pcr_n, "cbc", placeholder=True, inferred=inf), seg("s7", nx.S7, "s7")]),
    ]
    return "<h2>Key oligos</h2><seq>" + "".join(rows) + "</seq>"


def sci3_oligos() -> str:
    inf = False
    rows = [
        oligo("N5 ligation oligo", [seg("head", S.N5_HEAD, "r2", inferred=inf), seg("N5 barcode", "A" * 10, "cbc", placeholder=True, inferred=inf), seg("tail", S.N5_TAIL, "r2", inferred=inf)]),
        oligo("N7 ligation oligo", [seg("head", S.N7_HEAD, "r3", inferred=inf), seg("N7 barcode", "B" * 10, "cbc", placeholder=True, inferred=inf), seg("tail", S.N7_TAIL, "r3", inferred=inf)]),
        oligo("P5 PCR primer", [seg("P5", il.P5, "p5"), seg("i5", "I" * 10, "cbc", placeholder=True), seg("N5 head", S.N5_HEAD, "r2", inferred=inf)]),
        oligo("P7 PCR primer", [seg("P7", il.P7, "p7"), seg("i7", "I" * 10, "cbc", placeholder=True), seg("N7 head", S.N7_HEAD, "r3", inferred=inf)]),
    ]
    return "<h2>Key oligos</h2><seq>" + "".join(rows) + "</seq>"


def indexed_transposome(which: str) -> Scene:
    top = S.t5() if which == "T5" else S.t7()
    sc = Scene(); sc.strand(which, top, label=which); sc.anneal("pMENTS", [seg("ME'", nx.ME_RC, "me")], to=which, pair=("ME'", "ME"), label="pMENTS", mod5="p")
    sc.mark(which, "t5 barcode" if which == "T5" else "t7 barcode", "round-1 barcode")
    return sc


def two_tagmented() -> Scene:
    g = nx.TAGMENTATION_GAP
    top = [*S.t5(), seg("gap-facing A", "X" * g, placeholder=True), seg("core", "XXXXXXXX...XXXXXXXX", placeholder=True)]
    bot = [*S.t7(), seg("gap-facing B", "x" * g, placeholder=True), *complement_segments([top[-1]])]
    sc = Scene(); sc.strand("T5 end", top); sc.anneal("T7 end", bot, to="T5 end", pair=("core'", "core")); sc.anneal("T5 pMENTS", [seg("ME", nx.ME_RC, "me")], to="T5 end", pair=("ME", "ME"), label=""); sc.anneal("T7 pMENTS", [seg("ME", nx.ME_RC, "me")], to="T7 end", pair=("ME", "ME"), label="", above=True)
    sc.mark("T5 end", "gap-facing A", "9-nt gap"); sc.mark("T7 end", "gap-facing B", "9-nt gap")
    return sc


def final_panel(protocol: str) -> str:
    lib = S.library(protocol)
    caption = "Final four-part cell-barcode library."
    if protocol == S.SCI15: caption = "INFERRED — final library using the compatible 8-nt PCR-index layout; dotted indices are unavailable in the defining supplement."
    if protocol == S.SCI3: caption = "Final three-level library assembled from Supplementary Table S7."
    return panel([strand_row(lib), strand_row(lib, "bottom"), *annotation_rows(lib)], caption=caption)


def two_steps(protocol: str) -> str:
    return f'''<h2>Library generation</h2>
<h3>(1) Assemble indexed T5 and T7 transposomes</h3>{panel(indexed_transposome("T5").rows(), cls="small")}{panel(indexed_transposome("T7").rows(), cls="small")}
<h3>(2) Tagment nuclei in the first indexing plate</h3>{panel(two_tagmented().rows(), caption="Representative amplifiable T5...T7 product before gap fill; each end records the first-round well.")}
<h3>(3) Pool nuclei, redistribute, fill gaps and PCR-index</h3>{info('A P5/P7 primer pair records the second-round well. Same-end T5...T5 and T7...T7 fragments lack one primer site.')}
<h3>(4) Final library</h3>{final_panel(protocol)}'''


def ligation_scene(which: str) -> Scene:
    if which == "N5":
        joined = [seg("head", S.N5_HEAD, "r2"), seg("barcode", "A" * 10, "cbc", placeholder=True), seg("tail", S.N5_TAIL, "r2"), seg("s5 first 8", nx.S5[:8], "s5"), seg("s5 rest", nx.S5[8:], "s5"), seg("ME", nx.ME, "me")]
        bridge = [seg("splint", revcomp(S.N5_TAIL + nx.S5[:8]), "r2")]
        pair = ("splint", "tail")
    else:
        joined = [seg("head", S.N7_HEAD, "r3"), seg("barcode", "B" * 10, "cbc", placeholder=True), seg("tail", S.N7_TAIL, "r3"), seg("s7 first 8", nx.S7[:8], "s7"), seg("s7 rest", nx.S7[8:], "s7"), seg("ME", nx.ME, "me")]
        bridge = [seg("splint", revcomp(S.N7_TAIL + nx.S7[:8]), "r3")]
        pair = ("splint", "tail")
    sc = Scene(); sc.strand("joined strand", joined); sc.anneal("splint", bridge, to="joined strand", pair=pair, unpaired=())
    sc.mark("joined strand", "tail", "ligation nick", through=("s5 first 8" if which == "N5" else "s7 first 8"))
    return sc


def sci3_steps() -> str:
    return f'''<h2>Library generation</h2>
<h3>(1) Tagment nuclei with plain Nextera Tn5, then phosphorylate transferred 5' ends</h3>{info('T4 PNK creates the 5\'-phosphate required for ligation; no cell barcode is added by Tn5.')}
<h3>(2) First split-pool ligation: N5 barcode onto s5 ends</h3>{panel(ligation_scene("N5").rows(), cls="small", caption="An 8+8 splint bridges the N5 tail and s5.")}
<h3>(3) Second split-pool ligation: N7 barcode onto s7 ends</h3>{panel(ligation_scene("N7").rows(), cls="small", caption="The end-specific N7 splint analogously joins only s7 ends.")}
<h3>(4) Redistribute nuclei, fill gaps and add i5/i7 by PCR</h3>{info('Only an N5...N7 fragment carries both 15-nt PCR landing sites.')}
<h3>(5) Final library</h3>{final_panel(S.SCI3)}'''


def sequencing(protocol: str) -> str:
    lib = S.library(protocol)
    primers = S.SCI3_PRIMERS if protocol == S.SCI3 else S.TWO_LEVEL_PRIMERS
    if sp.verify(lib, primers): raise ValueError("read layout does not match final library")
    if protocol == S.SCI18:
        rows = [("Read 1", "50", "genomic insert"), ("Index 1", "8 + 27 dark + 10", "t7 barcode + PCR i7"), ("Index 2", "8 + 21 dark + 10", "t5 barcode + PCR i5"), ("Read 2", "50", "genomic insert")]
    elif protocol == S.SCI15:
        rows = [("Read 1", "not available", "genomic insert"), ("Index 1", "8 + 27 constant + PCR i7", "first- and second-round barcodes"), ("Index 2", "8 + 21 constant + PCR i5", "first- and second-round barcodes"), ("Read 2", "not available", "genomic insert")]
    else:
        rows = [("Read 1", "51", "genomic insert"), ("Index 1", "10 + 15 dark + 10", "N7 + i7"), ("Index 2", "10 + 15 dark + 10", "N5 + i5"), ("Read 2", "51", "genomic insert")]
    return ('<h2>Read layout</h2>' + sp.diagram(lib, primers)
            + table(("Read", "Cycles", "Content"), rows) + '</div>')


def render_page(protocol: str, out: Path) -> None:
    body = two_steps(protocol) if protocol != S.SCI3 else sci3_steps()
    olig = two_oligos(protocol) if protocol != S.SCI3 else sci3_oligos()
    out.write_text("\n".join([head(f"{protocol} library chemistry"), preamble(protocol), olig, body, sequencing(protocol)]), encoding="utf-8")
    print(f"wrote {out}  ({out.stat().st_size:,} bytes)")
