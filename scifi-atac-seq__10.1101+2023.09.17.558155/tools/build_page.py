#!/usr/bin/env python3
"""Build the diagram-centric scifi-ATAC-seq chemistry page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import nextera as nx
import scifi_atac as S
import seqprimers as sp
from chemdraw import annotation_rows, oligo, panel, protocol_graph_panels, strand_row
from page import head, info, table
from protocol_graph import AnnotationSlot
from protocol_views import assemble_page

OUT = HERE.parent / "scifi-atac-seq.html"


def oligos() -> str:
    rows = [
        oligo("Tn5 A adapter", S.adapter_a()),
        oligo("Tn5 B adapter", S.adapter_b()),
        oligo("Shared ME bottom", [S.seg("ME reverse complement", nx.ME_RC, "me")], mods="/5Phos/"),
        oligo("10x scATAC v1.1 gel-bead oligo", S.bead_oligo(), five="bead-|-5'-"),
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
{info('Four custom primers are used. Read 1 and Read 2 each begin with a 5-nt Tn5 barcode then the 19-nt ME; Index 1 reports sample index and Index 2 the 16-nt GEM barcode.')}
{sp.diagram(lib, S.SEQ_PRIMERS)}
{table(("Read", "Primer", "Primer site", "First bases"), rows)}
{sp.read_content_table(lib, S.SEQ_PRIMERS)}'''


def main() -> None:
    graph = S.protocol_graph()
    page = assemble_page(graph, [head("scifi-ATAC-seq library chemistry"), '<div class="wrap" data-protocol-complete="true">',
        '<h1>scifi-ATAC-seq &mdash; plate pre-indexing plus overloaded 10x scATAC</h1>',
        info('Defining source: <a href="https://doi.org/10.1101/2023.09.17.558155">Wang et al., bioRxiv (2023)</a>; droplet chemistry follows the cited 10x scATAC v1.1 guide.'),
        AnnotationSlot.INTRO, oligos(), construction(graph), sequencing(), '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
