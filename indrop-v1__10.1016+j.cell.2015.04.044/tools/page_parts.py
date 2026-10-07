#!/usr/bin/env python3
"""Shared renderer for individual inDrop v1 and v2 protocol pages."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import indrop as I
from chemdraw import (Construct, Scene, Segment, annotation_rows, complement_segments,
                      oligo, panel, revcomp, strand_row)
import illumina as il
from page import head, info, table
import seqprimers as sp

PAPERS = {
    I.V1: ("10.1016/j.cell.2015.04.044", "inDrop v1"),
    I.V2: ("10.1038/nprot.2016.154", "inDrop v2"),
}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def preamble(protocol: str) -> str:
    doi, title = PAPERS[protocol]
    caveat = ('''<div class="caveat"><b>Upstream-only v1 oligos.</b> The Cell supplement containing the author sequences was unavailable. Every dotted v1 oligo and final-library base follows the secondary scg_lib_structs reconstruction; only the reaction outline is established by the defining paper.</div>'''
              if protocol == I.V1 else "")
    return f'''<div class="wrap"><h1>{title} &mdash; droplet barcoding with photo-released hydrogel primers</h1>
{info(f'Defining source: <a href="https://doi.org/{doi}">doi:{doi}</a>.')}
{info('Each bead carries a split-pool cell barcode, a 6-nt UMI, oligo-dT and a T7 promoter. UV releases the primers into the droplet; pooled cDNA is amplified linearly by IVT.')}{caveat}'''


def bead_oligos() -> str:
    rows = [
        oligo("Acrydite bead primer", [seg("leader", I.LEADER), seg("T7 promoter", I.T7_PROMOTER, "tso"), seg("middle", I.MID), seg("PE1 partial", I.PE1[:24], "r2")], mods="Acrydite; photocleavable spacer"),
        oligo("P1 plate", [seg("W1*", I.W1_STAR, "r3"), seg("barcode 1", "A" * 8, "cbc", placeholder=True), seg("PE1*", revcomp(I.PE1), "r2")]),
        oligo("P2 plate", [seg("BA19", "B" + "A" * 19, placeholder=True), seg("N6", "N" * 6, "umi", placeholder=True), seg("barcode 2", "B" * 8, "cbc", placeholder=True), seg("W1*", I.W1_STAR, "r3")]),
    ]
    return "<h2>Bead synthesis oligos</h2><seq>" + "".join(rows) + "</seq>"


def library_oligos(protocol: str) -> str:
    inf = protocol == I.V1
    if protocol == I.V1:
        rows = [
            oligo("RNA ligation oligo (RLO)", [seg("RLO", I.RLO, "r2", inferred=inf)], mods="5Phos; 3SpC3"),
            oligo("Second RT primer", [seg("RT2", "GTCTCGGCATTCCTGCTGAAC", "r2", inferred=inf)]),
            oligo("PCR enrichment primer 1", [seg("P5", il.P5, "p5", inferred=inf), seg("PE1", "TCTTTCCCTACACGA", "r1", inferred=inf)]),
            oligo("PCR enrichment primer 2", [seg("P7", il.P7, "p7", inferred=inf), seg("RT2", "CGGTCTCGGCATTCCTGCTGAAC", "r2", inferred=inf)]),
        ]
    else:
        rows = [
            oligo("PE2-N6", [seg("PE2", I.PE2, "r1"), seg("N6", "N" * 6, "umi", placeholder=True)]),
            oligo("PE2 PCR primer", [seg("P5", il.P5, "p5"), seg("GGTC", "GGTC"), seg("PE2 anneal", I.PE2[:18], "r1")]),
            oligo("PE1 indexed PCR primer", [seg("P7", il.P7, "p7"), seg("i7 as ordered", "I" * 6, "cbc", placeholder=True), seg("PE1 anneal", I.PE1[:16], "r2")]),
            oligo("Custom Read 1 / Index / Read 2", [seg("R1", I.V2_R1, "r1"), seg(" / ", "", placeholder=True), seg("Index", I.V2_INDEX, "r2"), seg(" / ", "", placeholder=True), seg("R2", I.V2_R2, "r2")]),
        ]
    return "<h2>Library oligos</h2><seq>" + "".join(rows) + "</seq>"


def bead_steps(protocol: str) -> str:
    bead = I.bead_primer(inferred=(protocol == I.V1))
    caption = "Finished bead oligo; barcode 1 is 8-11 nt, so one representative length is drawn."
    if protocol == I.V1:
        caption = "INFERRED — finished v1 bead oligo. Dotted bases use the v2 supplementary sequences and secondary v1 reconstruction; barcode 1 is 8-11 nt."
    return f'''<h2>Bead and droplet steps</h2>
<h3>(1) Build the bead barcode in two 384-way split-pool extensions</h3>{panel([strand_row(bead), *annotation_rows(bead)], caption=caption)}
<h3>(2) Co-encapsulate one cell and one bead; release primers with UV</h3>{info('The photocleavable spacer frees many copies of one bead barcode into the droplet.')}
<h3>(3) Reverse-transcribe from T19V, then make the T7 promoter double-stranded</h3>{panel([strand_row(Construct([*bead.segments, seg("antisense cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True, inferred=(protocol == I.V1))]))], caption=("INFERRED — v1 first-strand structure using the reconstructed bead oligo." if protocol == I.V1 else "First strand carries bead barcode and UMI upstream of antisense cDNA."))}
<h3>(4) T7 in-vitro transcription</h3>{info('IVT linearly amplifies the pooled material as antisense RNA carrying the UMI and both bead-barcode blocks.')}
'''


def v1_steps() -> str:
    lig = Construct([seg("aRNA fragment", "XXXXXXXX...XXXXXXXX", placeholder=True, inferred=True), seg("RLO", I.RLO, "r2", inferred=True)])
    return f'''<h2>v1 library conversion</h2>
<h3>(5) Fragment aRNA and ligate the blocked RLO</h3>{panel([strand_row(lig), *annotation_rows(lig)], caption="INFERRED — secondary-source v1 RLO ligated to the aRNA 3' end.")}
<h3>(6) Second RT from RLO, then P5/P7 enrichment PCR</h3>{final_panel(I.V1)}'''


def v2_steps() -> str:
    rt = Construct([seg("PE2", I.PE2, "r1"), seg("N6", "N" * 6, "umi", placeholder=True), seg("sense cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True), seg("bead-meta copy", "AAAA...AAAA", placeholder=True)])
    return f'''<h2>v2 library conversion</h2>
<h3>(5) Random-prime the aRNA with PE2-N6</h3>{panel([strand_row(rt), *annotation_rows(rt)], caption="PE2-N6 makes sense-oriented cDNA from the antisense aRNA; the opposite end retains the bead metadata.")}
<h3>(6) Add P5 on the PE2 end and indexed P7 on the PE1 end</h3>{final_panel(I.V2)}'''


def final_panel(protocol: str) -> str:
    lib = I.library(protocol)
    caption = "Final library; one representative 8-nt barcode-1 length is drawn."
    if protocol == I.V1:
        caption = "INFERRED — entire v1 final construct is dotted because its oligos are available only from the secondary source."
    return panel([strand_row(lib), strand_row(lib, "bottom"), *annotation_rows(lib)], caption=caption)


def sequencing(protocol: str) -> str:
    lib = I.library(protocol); primers = I.V1_PRIMERS if protocol == I.V1 else I.V2_PRIMERS
    if any(sp.locate(lib, p) is None for p in primers): raise ValueError("inDrop read site missing")
    rows = ([("Read 1", "51 (secondary source)", "barcode 1 &middot; W1 &middot; barcode 2 &middot; UMI &middot; poly(T)"), ("Read 2", "secondary source", "cDNA")]
            if protocol == I.V1 else
            [("Read 1", "not available", "sense cDNA"), ("Index 1", "6", "sample index"), ("Read 2", "at least 51", "barcode 1 &middot; W1 &middot; barcode 2 &middot; UMI &middot; poly(T)")])
    return '<h2>Read layout</h2>' + table(("Read", "Cycles", "Content"), rows) + '</div>'


def render_page(protocol: str, out: Path) -> None:
    version_steps = v1_steps() if protocol == I.V1 else v2_steps()
    out.write_text("\n".join([head(f"{protocol} library chemistry"), preamble(protocol), bead_oligos(), library_oligos(protocol), bead_steps(protocol), version_steps, sequencing(protocol)]), encoding="utf-8")
    print(f"wrote {out}  ({out.stat().st_size:,} bytes)")
