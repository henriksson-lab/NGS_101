#!/usr/bin/env python3
"""Build the diagram-centric HyDrop-RNA chemistry page."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import hydrop_rna as H
import seqprimers as sp
from chemdraw import annotation_rows, oligo, panel, protocol_graph_panels, strand_row
from page import head, info

OUT = HERE.parent / "hydrop-rna.html"


def oligos() -> str:
    rows = [
        oligo("Final split-pool-barcoded RNA bead oligo", list(H.bead_oligo()),
              mods="/5Acryd//iThioMC6-D/"),
        oligo("Template-switching oligo", H.tso_oligo()),
        oligo("TSO-P cDNA-PCR primer", [H.seg("SMART handle", H.TSO_P, "tso")]),
        oligo("HYi5_TruSeq library-PCR primer", H.hyi5_primer()),
        oligo("HYi7 library-PCR primer", H.hyi7_primer(),
              mods="last two bonds are phosphorothioate"),
    ]
    return "<h2>Key oligos</h2><seq>" + "".join(rows) + "</seq>"


def construction() -> str:
    return '<h2>Library construction</h2>\n'+protocol_graph_panels(H.protocol_graph(),cls="long")


def sequencing() -> str:
    lib = H.final_library()
    return sp.section(lib, H.SEQ_PRIMERS, heading="Sequencing",
        intro='Read 1: 50 cycles of cDNA; Index 1 and Index 2: 10-cycle sample indices; Read 2: 58 cycles covering the three cell barcodes, two linkers and 8-nt UMI.',
        read_lengths={"Read 1":50,"Index 1 (i7)":10,"Index 2 (i5)":10,"Read 2":58},
        read_uses={"Read 2":"extract"}, show_read_content=True)


def main() -> None:
    page = "\n".join([head("HyDrop-RNA library chemistry"), '<div class="wrap" data-protocol-complete="true">',
        '<h1>HyDrop-RNA &mdash; droplet single-cell RNA sequencing</h1>',
        info('Defining source: <a href="https://doi.org/10.7554/eLife.73971">De Rop et al., <i>eLife</i> (2022)</a>.'),
        oligos(), construction(), sequencing(), '</div>'])
    OUT.write_text(page, encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
