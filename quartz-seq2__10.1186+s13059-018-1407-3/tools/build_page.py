#!/usr/bin/env python3
"""Build the diagram-centric Quartz-Seq2 chemistry page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import quartz_seq2 as Q
import seqprimers as sp
from chemdraw import Construct, Segment, oligo, protocol_graph_panels
from page import head, info, table


OUT = HERE.parent / "quartz-seq2.html"


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


def oligos() -> str:
    rows = [
        oligo("1536-well RT primer", Q.rt_primer()),
        oligo("Tagging primer", Q.tagging_primer()),
        oligo("gM suppression primer", [seg("G", "G"), seg("M handle", Q.M, "tso")]),
        oligo("rYshapeP5 short strand", [seg("stem", Q.RYSHAPE_P5[:12], "r1"),
                                          seg("unpaired arm", Q.RYSHAPE_P5[12:], "r1")]),
        oligo("rYshapeP7LT06 long strand",
              [seg("P7", Q.RYSHAPE_P7_LT06[:24], "p7"),
               seg("i7", Q.LT06_OLIGO_INDEX, "cbc"),
               seg("TruSeq Read 2", Q.RYSHAPE_P7_LT06[-34:], "r2")]),
        oligo("P5-gMac selective PCR primer",
              [seg("truncated P5", Q.P5_GMAC[:28], "p5"), seg("TT", "TT"),
               seg("gM + AC", Q.P5_GMAC[30:], "tso")],
              mods="two 3'-terminal phosphorothioates"),
        oligo("Read1DropQuartz sequencing primer", [seg("custom R1", Q.READ1_DROPQUARTZ, "r1")]),
    ]
    return "<h2>Key oligos</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>"



def sequencing() -> str:
    lib = Q.final_library()
    rows = []
    for primer, hit in Q.primer_landings():
        rows.append((html.escape(primer.role), html.escape(primer.name),
                     html.escape(", ".join(hit.covers)),
                     f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('Quartz-Seq2 v3.2 uses Read 1: 23 cycles, Index 1: 6, and Read 2: 63. Read 1 reports the 15-base cell barcode followed by the 8-base UMI; Read 2 reports cDNA.')}
{sp.diagram(lib, Q.SEQ_PRIMERS)}
{table(("Read", "Primer", "Primer site", "First bases"), rows)}
<h3>Read 1 barcode layout</h3>
{table(("Cycles", "Content"), Q.read1_layout(), scroll=False)}
{sp.read_content_table(lib, Q.SEQ_PRIMERS,
    {"Read 1":23,"Index 1 (i7)":6,"Read 2":63},
    acquired_roles=tuple(p.role for p in Q.SEQ_PRIMERS))}
'''


def main() -> None:
    graph = Q.protocol_graph()
    construction = "<h2>Library construction</h2>\n" + protocol_graph_panels(graph)
    page = "\n".join([head("Quartz-Seq2 library chemistry"),
        '<div class="wrap" data-protocol-complete="true">',
        '<h1>Quartz-Seq2 &mdash; barcoded 3-prime single-cell RNA sequencing</h1>',
        info('Defining source: <a href="https://doi.org/10.1186/s13059-018-1407-3">Sasagawa et al., <i>Genome Biology</i> (2018)</a>.'),
        info('Published-source discrepancy: Additional file 5, Table S4, v3.2 sheet '
             '(MOESM5 extracted lines 389–393) labels its heading “SeqLv15mer” but the '
             'next line calls the same v3.2 barcodes “14 mer.” The printed eMDRT0001 '
             'barcode <code>ATCAATCCTATCTGC</code> has 15 bases, as do the v3.2 '
             'barcodes inspected in that sheet; the article Methods (cell-barcode design) '
             'also specifies 15 nt for v3.2. This model uses 15 nt. The conflicting '
             '“14 mer” line appears to be a table heading error; no correction notice '
             'was found in the inspected sources.'),
        info('RNase H and TdT act in the same source tailing mixture. Native second-strand and adapter-ligated precursors enter the two PCR stages; residual original polymers remain explicitly bounded.'), oligos(), construction, sequencing(), '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
