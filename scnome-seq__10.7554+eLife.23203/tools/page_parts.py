#!/usr/bin/env python3
"""Shared renderer for individual scNOMe-seq and scCOOL-seq pages."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import nomecool as N
from chemdraw import Construct, Scene, Segment, annotation_rows, complement_segments, oligo, panel, strand_row
import illumina as il
from page import head, info, table
import seqprimers as sp

PAPERS = {
    N.SCNOME: ("10.7554/eLife.23203", "scNOMe-seq"),
    N.SCCOOL: ("10.1038/cr.2017.82", "scCOOL-seq"),
}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def preamble(protocol: str) -> str:
    doi, title = PAPERS[protocol]
    notes = ("01_scnome-seq-sccool-seq.html" if protocol == N.SCNOME else
             "../scnome-seq__10.7554+eLife.23203/01_scnome-seq-sccool-seq.html")
    if protocol == N.SCNOME:
        caveat = '<div class="caveat"><b>Commercial kit interior.</b> Pico Methyl-Seq random primers and internal adapters are proprietary. Dotted question marks preserve the unknown insert junctions; the published P5/P7 amplification primers fix the outer library arms.</div>'
    else:
        caveat = '<div class="caveat"><b>PCR-primer sequences unavailable.</b> The paper prints both random primers but only names Forward PE1.0 and the NEB indexed reverse primer. Dotted outer arms show the standard TruSeq-shaped model implied by those reagents.</div>'
    return f'''<div class="wrap"><h1>{title} &mdash; accessibility and DNA methylation from one cell</h1>
<p class="research-notes"><a href="{notes}">Research notes</a></p>
{info(f'Defining source: <a href="https://doi.org/{doi}">doi:{doi}</a>.')}
{info('M.CviPI marks accessible GpC cytosines before bisulfite conversion. Retained GpC-C reports accessibility; retained CpG-C reports endogenous DNA methylation.')}{caveat}'''


def oligos(protocol: str) -> str:
    if protocol == N.SCNOME:
        rows = [
            oligo("forward_p5", [seg("P5", il.P5, "p5"), seg("Read 1 remainder", il.TRUSEQ_READ1[4:], "r1")]),
            oligo("INDEXnn", [seg("P7", il.P7, "p7"), seg("6-nt index", "I" * 6, "cbc", placeholder=True), seg("Read 2 arm", il.TRUSEQ_READ2, "r2")]),
            oligo("Pico Methyl-Seq internal primers", [seg("sequence unavailable", "?", inferred=True, placeholder=True)]),
        ]
    else:
        rows = [
            oligo("Random primer 1", N.primer1_segments(), mods="5'-biotin"),
            oligo("Random primer 2", N.primer2_segments()),
            oligo("Forward PE1.0 (sequence inferred from reagent name)", [seg("P5", il.P5, "p5", inferred=True), seg("Read 1 remainder", il.TRUSEQ_READ1[4:], "r1", inferred=True)]),
            oligo("NEB reverse indexed primer", [seg("sequence unavailable", "?", inferred=True, placeholder=True)]),
        ]
    return "<h2>Key oligos</h2><seq>" + "".join(rows) + "</seq>"


def mark_and_convert() -> str:
    before = Construct([seg("CpG", "CG", placeholder=True), seg("spacer", "...", placeholder=True), seg("open GpC", "GC", placeholder=True), seg("spacer2", "...", placeholder=True), seg("protected GpC", "GC", placeholder=True)])
    marked = Construct([seg("endogenous mCpG", "CG", placeholder=True), seg("spacer", "...", placeholder=True), seg("M.CviPI-marked GpC", "GC", placeholder=True), seg("spacer2", "...", placeholder=True), seg("unmarked GpC", "GC", placeholder=True)])
    converted = Construct([seg("retained CpG-C", "CG", placeholder=True), seg("spacer", "...", placeholder=True), seg("retained open GpC-C", "GC", placeholder=True), seg("spacer2", "...", placeholder=True), seg("converted protected GpT", "GT", placeholder=True)])
    return f'''<h2>Molecular marking and conversion</h2>
<h3>(1) Mark accessible GpC with M.CviPI</h3>{panel([strand_row(before), *annotation_rows(before), strand_row(marked), *annotation_rows(marked)], cls="small", caption="M.CviPI methylates exposed GpC; nucleosome- or factor-protected GpC remains unmethylated.")}
<h3>(2) Bisulfite-convert the DNA</h3>{panel([strand_row(converted), *annotation_rows(converted)], cls="small", caption="Unmethylated cytosine becomes U/T in the library; methylated CpG and enzyme-marked GpC remain C.")}'''


def cool_priming() -> str:
    first = Construct([*N.primer1_segments(), seg("copied converted DNA", "XXXXXXXX...XXXXXXXX", placeholder=True)])
    top = [seg("Read 1 handle", N.COOL_P1_HANDLE, "r1"), seg("primer-1 N9", "N" * 9, placeholder=True), seg("insert", "XXXXXXXX...XXXXXXXX", placeholder=True), seg("primer-2 site", "X" * 9, placeholder=True)]
    bottom = [*N.primer2_segments(), *complement_segments(top[:-1])]
    sc = Scene(); sc.strand("captured first strand", top, mod5="biotin"); sc.anneal("second strand", bottom, to="captured first strand", pair=("Read 1 handle'", "Read 1 handle"))
    return f'''<h2>Library generation</h2>
<h3>(3) Random-prime the first strand</h3>{panel([strand_row(first), *annotation_rows(first)], caption="Biotinylated primer 1 adds the Read-1-side handle and an N9 priming tract.")}
<h3>(4) Capture, remove the converted template, and random-prime the second strand</h3>{panel(sc.rows(), caption="Primer 2 adds a different handle and N9 on the bead-bound first-strand product.")}'''


def nome_library_steps() -> str:
    return f'''<h2>Library generation</h2>
<h3>(3) Build the bisulfite library with Pico Methyl-Seq</h3>{info('The low-input kit random-primes converted single-stranded DNA and pre-amplifies it. Its internal oligo sequences are not published.')}
<h3>(4) Add the published P5 and 6-nt indexed P7 primers</h3>{final_panel(N.SCNOME)}'''


def final_panel(protocol: str) -> str:
    lib = N.library(protocol)
    reason = ("Pico Methyl-Seq insert junctions" if protocol == N.SCNOME else
              "the two unprinted PCR primers")
    return panel([strand_row(lib), strand_row(lib, "bottom"), *annotation_rows(lib)],
                 caption=f"INFERRED — final paired-end library. Dotted bases represent {reason}.")


def sequencing(protocol: str) -> str:
    lib = N.library(protocol)
    if any(sp.locate(lib, p) is None for p in N.SEQ_PRIMERS):
        raise ValueError("sequencing primer missing")
    rows = ([("Read 1", "100", "bisulfite insert; first 6 bases clipped"), ("Index 1", "6", "cell/sample index"), ("Read 2", "100", "bisulfite insert; non-directional")]
            if protocol == N.SCNOME else
            [("Read 1", "150", "primer-1 N9, then bisulfite insert"), ("Index 1", "indexed", "cell/sample index"), ("Read 2", "150", "primer-2 N9, then bisulfite insert")])
    return ('<h2>Read layout</h2>' + sp.diagram(lib, N.SEQ_PRIMERS)
            + table(("Read", "Cycles", "Content"), rows)
            + info('There is no molecular cell barcode or UMI; the final sample index identifies the well.')
            + '</div>')


def render_page(protocol: str, out: Path) -> None:
    lib_steps = nome_library_steps() if protocol == N.SCNOME else cool_priming() + '<h3>(5) Indexed PCR on the beads</h3>' + final_panel(N.SCCOOL)
    out.write_text("\n".join([head(f"{protocol} library chemistry"), preamble(protocol), oligos(protocol), mark_and_convert(), lib_steps, sequencing(protocol)]), encoding="utf-8")
    print(f"wrote {out}  ({out.stat().st_size:,} bytes)")
