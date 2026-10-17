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
from protocol_graph import AnnotationSlot
from protocol_views import assemble_page

OUT = HERE.parent / "scifi-rna-seq.html"


def oligos() -> str:
    rows = [
        oligo("Round-1 barcoded RT primer", S.rt_primer(), mods="/5Phos/"),
        oligo("10x scATAC gel-bead oligo", S.bead_oligo(), five="bead-|-5'-"),
        oligo("Bridge oligo", S.bridge_oligo(), three="-[ddC]-3'"),
        oligo("Template-switching oligo", S.tso_oligo()),
        oligo("TSO_enrichment_primer", [S.seg("SMART handle", S.rt.SMART_HANDLE, "tso")]),
        oligo("Partial-P5", [S.seg("partial P5", S.PARTIAL_P5, "p5")]),
        oligo("i7-only Tn5 adaptor (Tn5-bottom_Read2N)", [S.seg("s7", nx.S7, "s7"),
                                                          S.seg("mosaic end", nx.ME, "me")]),
        oligo("Tn5-top_ME", [S.seg("mosaic end reverse complement", nx.ME_RC, "me")],
              mods="/5Phos/"),
        oligo("Indexed P7 primer (P7-001 to P7-016)", list(S.indexed_p7_primer())),
    ]
    return "<h2>Key oligos</h2><seq>" + "".join(rows) + "</seq>"


def construction(graph) -> str:
    return '<h2>Library construction</h2>\n'+protocol_graph_panels(graph,cls="long")


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
    graph = S.protocol_graph()
    page = assemble_page(graph, [head("scifi-RNA-seq library chemistry"), '<div class="wrap" data-protocol-complete="true">',
        '<h1>scifi-RNA-seq &mdash; combinatorial fluidic indexing</h1>',
        info('Defining source: <a href="https://doi.org/10.1038/s41592-021-01153-z">Datlinger et al., Nature Methods (2021)</a>, its Supplementary Table 1 and Supplementary Protocol v2021-01-25. The gel-bead oligo is printed in no table; it is drawn only in the published Extended Data Fig. 3a (an image) and read from that rendered figure.'),
        AnnotationSlot.INTRO, oligos(), construction(graph), sequencing(), '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
