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
from protocol_graph import AnnotationSlot
from protocol_views import assemble_page

OUT = HERE.parent / "sci-rna-seq.html"


def oligos() -> str:
    rows = [oligo("Barcoded anchored oligo-dT RT primer", S.rt_primer()),
            oligo("Indexed P5 PCR primer", S.p5_primer()),
            oligo("Indexed P7 PCR primer", S.p7_primer())]
    return "<h2>Key oligos</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>"


def construction(graph) -> str:
    return '<h2>Library construction</h2>\n'+protocol_graph_panels(graph,cls="long")


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
    graph = S.protocol_graph()
    page=assemble_page(graph, [head("sci-RNA-seq library chemistry"),
        '<div class="wrap" data-protocol-complete="true">',
        '<h1>sci-RNA-seq &mdash; two-level combinatorial indexing of fixed cells</h1>',
        info('Defining source: <a href="https://doi.org/10.1126/science.aam8940">Cao et al., <i>Science</i> (2017)</a>. '
             'The displayed implementation follows the original <a href="https://doi.org/10.1101/104844">2017 preprint Methods</a>, '
             'which print the RT/PCR primer templates and NextSeq 18/10/10/52 run. The final Science supplement has not yet been recovered, '
             'so final-version comparison remains an open source-review task. The original <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98561">GEO annotations</a> supply experiment-specific RT barcode lists; the displayed index placeholders leave individual well assignments unselected. '
             'The RT barcode and PCR-well i5/i7 pair jointly identify a cell; barcode collisions remain possible. '
             'INFERRED — the selected gap-filled intermediate uses a canonical Nextera mosaic end; exact TDE1 intermediate structures are unassigned. '
             'This page models the two-level branch. The final article also tests a third level with indexed Tn5; '
             'that branch differs from the later 2019 sci-RNA-seq3 ligation workflow.'),
        AnnotationSlot.INTRO, oligos(), construction(graph), AnnotationSlot.AFTER_WORKFLOW,
        sequencing(), '</div>'])
    OUT.write_text(page,encoding="utf-8"); print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
