#!/usr/bin/env python3
"""Build the diagram-centric Microwell-seq chemistry page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import microwell_seq as M
import nextera as nx
import seqprimers as sp
from chemdraw import annotation_rows, oligo, panel, strand_row
from page import head, info, table

OUT = HERE.parent / "microwell-seq.html"


def oligos() -> str:
    rows = [
        oligo("Final seqA split-pool bead oligo", list(M.bead_oligo()), five="bead-|-5'-"),
        oligo("Template-switching oligo", M.tso_oligo()),
        oligo("TSO-PCR primer", [M.seg("SMART handle", M.rt.SMART_HANDLE, "tso")]),
        oligo("P5 bead-end-selective primer", [M.seg("P5", M.il.P5, "p5"),
                                                M.seg("spacer", M.P5_SPACER),
                                                M.seg("SMART handle", M.rt.SMART_HANDLE, "tso"),
                                                M.seg("AC", "AC")], three="-*A*C-3'"),
        oligo("N701 library-PCR primer", [M.seg("P7", M.il.P7, "p7"),
                                           M.seg("i7", M.I7_OLIGO_INDEX, "cbc"),
                                           M.seg("s7", nx.S7, "s7")]),
    ]
    return "<h2>Key oligos</h2><seq>" + "".join(rows) + "</seq>"


def construction() -> str:
    bead = M.bead_oligo()
    lib = M.final_library()
    return f'''<h2>Library construction</h2>
<h3>(1) Build the bead barcode by three split-pool extension rounds</h3>
{panel([strand_row(bead), *annotation_rows(bead)], cls="small",
       caption="The seqA bead carries three 6-nt barcode blocks, two fixed 15-nt linkers, a 6-nt UMI and T30.")}
<h3>(2) Capture microwell RNA, reverse-transcribe and template-switch</h3>
{panel(M.rt_scene().rows(), cls="small",
       caption="Collected beads are pooled for reverse transcription; the TSO installs the second SMART handle.")}
<h3>(3) Single-primer WTA</h3>
{panel([strand_row(bead), *annotation_rows(bead)], cls="small",
       caption="TSO-PCR amplifies SMART-handle-flanked cDNA before library construction.")}
<h3>(4) INFERRED &mdash; symmetric s7 tagmentation and bead-end-selective PCR</h3>
{panel([*M.final_scene().rows(), *annotation_rows(lib)], cls="small",
       caption="INFERRED — the defining paper prints the P5 and indexed-P7 primers but not the Tn5 adaptor. Dotted s7/ME and P7-side regions follow those primers and the later same-lab description.")}'''


def sequencing() -> str:
    lib = M.final_library()
    rows = []
    for primer in M.SEQ_PRIMERS:
        hit = sp.locate(lib, primer)
        rows.append((html.escape(primer.role), html.escape(primer.name),
                     html.escape(", ".join(hit.covers)),
                     f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('Read 1 needs 54 cycles to cover the three cell-barcode blocks and UMI; Index 1 reports the sample index; Read 2 reports cDNA.')}
{sp.diagram(lib, M.SEQ_PRIMERS)}
{table(("Read", "Primer", "Primer site", "First bases"), rows)}
<h3>Read 1 layout</h3>
{table(("Cycles", "Content"), M.read1_layout(), scroll=False)}'''


def main() -> None:
    page = "\n".join([head("Microwell-seq library chemistry"), '<div class="wrap">',
        '<h1>Microwell-seq &mdash; split-pool-barcoded bead capture</h1>',
        info('Defining source: <a href="https://doi.org/10.1016/j.cell.2018.02.001">Han et al., <i>Cell</i> (2018)</a>.'),
        '<div class="caveat"><b>Library-adaptor boundary.</b> The defining source does not print the Tn5 adaptor; inferred regions are dotted.</div>',
        oligos(), construction(), sequencing(), '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
