#!/usr/bin/env python3
"""Build the diagram-centric CH-ATAC-seq page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import ch_atac as C
import seqprimers as sp
from chemdraw import Scene, annotation_rows, oligo, panel
from page import head, info, table

OUT = HERE.parent / "ch-atac-seq.html"


def render() -> str:
    lib = C.final_library()
    duplex = Scene.duplex(list(lib), label="library")
    landings = [(html.escape(p.role), html.escape(p.name),
                 html.escape(", ".join(hit.covers)),
                 f"<code>{html.escape(hit.reads)}</code>&hellip;")
                for p, hit in C.primer_landings(lib)]
    oligos = "\n".join([
        oligo("Indexed Tn5 transfer strand", C.indexed_tn5_transfer()),
        oligo("Tn5 Primer C", C.primer_c()),
        oligo("HY head", C.hy_head()),
        oligo("Barcoded HY oligo", C.barcoded_hy()),
    ])
    return "\n".join([
        head("CH-ATAC-seq library chemistry"), '<div class="wrap">',
        '<h1>CH-ATAC-seq</h1>',
        info('Defining source: <a href="https://doi.org/10.1016/j.devcel.2024.01.015">Zhang et al., <i>Developmental Cell</i> (2024)</a>.'),
        '<h2>Three indexing rounds</h2>',
        info('Tn5 barcode (384 wells) &rarr; hybridization barcode (768 wells) &rarr; indexed PCR (96 wells). Pools are split again between rounds.'),
        '<h2>Key oligos</h2>', '<seq>', oligos, '</seq>',
        '<h2>Library construction</h2>',
        '<h3>(1) Tagment with indexed-handle and Primer-C transposomes</h3>',
        info('Productive fragments receive the SMART handle and Tn5 barcode at one end, and s7 from Primer C at the other.'),
        '<h3>(2) Hybridize the second barcode</h3>',
        info('A pre-annealed HY duplex uses its SMART-handle overhang to capture the indexed Tn5 end. The HY head supplies s5; a blocking oligo occupies free capture overhangs.'),
        '<h3>(3) Fill, amplify and add MGI ends</h3>',
        panel([*duplex.rows(), *annotation_rows(lib)], cls="small", caption="Final DNBSEQ library assembled from the defining paper's Supplementary Table 1."),
        '<h2>Sequencing</h2>',
        info('Paired 100-cycle reads and a 10-cycle index read are used, with Read 1 dark cycles over the constant handle and mosaic end.'),
        sp.diagram(lib, C.SEQ_PRIMERS),
        table(("Read", "Primer", "Primer site", "First bases"), landings),
        '<h3>Read 1 layout</h3>', table(("Cycles", "Content"), C.read1_layout(), scroll=False),
        '</div>'
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
