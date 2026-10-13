#!/usr/bin/env python3
"""Build the diagram-centric scifi-RNA-seq chemistry page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import nextera as nx
import scifi_rna as S
import seqprimers as sp
from chemdraw import annotation_rows, oligo, panel, protocol_graph_panels, strand_row
from page import head, info, table

OUT = HERE.parent / "scifi-rna-seq.html"


def oligos() -> str:
    rows = [
        oligo("Round-1 barcoded RT primer", S.rt_primer(), mods="/5Phos/"),
        oligo("10x scATAC gel-bead oligo", S.bead_oligo(), five="bead-|-5'-"),
        oligo("Bridge oligo", S.bridge_oligo(), three="-[ddC]-3'"),
        oligo("Template-switching oligo", [S.seg("SMART handle + GAAT", S.TSO[:-3], "tso"),
                                             S.seg("rGrGrG", "rGrGrG", placeholder=True)]),
        oligo("Partial-P5", [S.seg("partial P5", S.PARTIAL_P5, "p5")]),
        oligo("i7-only Tn5 adaptor", [S.seg("s7", nx.S7, "s7"),
                                       S.seg("mosaic end", nx.ME, "me")]),
    ]
    return "<h2>Key oligos</h2><seq>" + "".join(rows) + "</seq>"


def construction() -> str:
    return '<h2>Library construction</h2>\n'+protocol_graph_panels(S.protocol_graph(),cls="long")


def sequencing() -> str:
    lib = S.final_library()
    rows = []
    for primer in S.SEQ_PRIMERS:
        hit = sp.locate(lib, primer)
        rows.append((html.escape(primer.role), html.escape(primer.name),
                     html.escape(", ".join(hit.covers)),
                     f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('NovaSeq layout: Read 1 21 cycles; Index 1 8; Index 2 16; Read 2 78. Read 1 carries the plate barcode and UMI; Index 2 carries the bead barcode.')}
{sp.diagram(lib, S.SEQ_PRIMERS)}
{table(("Read", "Primer", "Primer site", "First bases"), rows)}
<h3>Read 1 layout</h3>
{table(("Cycles", "Content"), S.read1_layout(), scroll=False)}
{sp.read_content_table(lib,S.SEQ_PRIMERS,
    {"Read 1":21,"Index 1 (i7)":8,"Index 2 (i5)":16,"Read 2":78})}'''


def main() -> None:
    page = "\n".join([head("scifi-RNA-seq library chemistry"), '<div class="wrap" data-protocol-complete="true">',
        '<h1>scifi-RNA-seq &mdash; combinatorial fluidic indexing</h1>',
        info('Defining source: <a href="https://doi.org/10.1101/2019.12.17.879304">Datlinger et al., bioRxiv (2019)</a>; expanded protocol: <a href="https://doi.org/10.1038/s41592-021-01153-z">Nature Methods (2021)</a>.'),
        oligos(), construction(), sequencing(), '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
