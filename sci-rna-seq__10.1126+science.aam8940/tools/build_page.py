#!/usr/bin/env python3
"""Build the diagram-centric sci-RNA-seq chemistry page."""
from __future__ import annotations
import html
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "lib"))
import sci_rna_seq as S
import seqprimers as sp
from chemdraw import oligo, protocol_graph_panels
from page import head, info, table

OUT = HERE.parent / "sci-rna-seq.html"


def oligos() -> str:
    rows = [oligo("Barcoded anchored oligo-dT RT primer", S.rt_primer()),
            oligo("Indexed P5 PCR primer", S.p5_primer()),
            oligo("Indexed P7 PCR primer", S.p7_primer())]
    return "<h2>Key oligos</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>"


def construction() -> str:
    return '<h2>Library construction</h2>\n'+protocol_graph_panels(S.protocol_graph(),cls="long")


def sequencing() -> str:
    lib=S.final_library()
    rows=[]
    for p,h in S.primer_landings():
        rows.append((html.escape(p.role),html.escape(p.name),html.escape(", ".join(h.covers)),f"<code>{html.escape(h.reads)}</code>&hellip;"))
    return f'''<h2>Sequencing</h2>
{info('Read 1: 18 cycles; Index 1: 10; Index 2: 10; Read 2: 52. Read 1 contains UMI + RT barcode; Read 2 contains transcript sequence.')}
{sp.diagram(lib,S.SEQ_PRIMERS)}
{table(("Read","Primer","Primer site","First bases"),rows)}
<h3>Read 1 barcode layout</h3>{table(("Cycles","Content"),S.read1_layout(),scroll=False)}
{sp.read_content_table(lib, S.SEQ_PRIMERS,
    {"Read 1":18,"Index 1 (i7)":10,"Index 2 (i5)":10,"Read 2":52})}
'''


def main() -> None:
    page="\n".join([head("sci-RNA-seq library chemistry"),'<div class="wrap" data-protocol-complete="true">','<h1>sci-RNA-seq &mdash; two-level combinatorial indexing of fixed cells</h1>',info('Defining source: <a href="https://doi.org/10.1126/science.aam8940">Cao et al., <i>Science</i> (2017)</a>.'),oligos(),construction(),sequencing(),'</div>'])
    OUT.write_text(page,encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
