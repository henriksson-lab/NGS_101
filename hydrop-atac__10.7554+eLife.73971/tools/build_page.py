#!/usr/bin/env python3
"""Build the diagram-centric HyDrop-ATAC chemistry page."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import hydrop_atac as H
import seqprimers as sp
from chemdraw import Scene, annotation_rows, oligo, panel, protocol_graph_panels, strand_row
from page import head, info

OUT = HERE.parent / "hydrop-atac.html"


def oligos() -> str:
    bead = H.bead_oligo()
    rows = [
        oligo("Final split-pool-barcoded bead oligo", list(bead),
              mods="/5Acryd//iThioMC6-D/"),
        oligo("HYi7 bulk-PCR primer", H.hyi7_primer(), mods="last two bonds are phosphorothioate"),
        oligo("HYi5 bulk-PCR primer", H.hyi5_primer()),
    ]
    return ("<h2>Bead oligo and PCR primers</h2>\n<seq>\n" + "\n".join(rows) + "\n</seq>\n"
            + panel([strand_row(bead), *annotation_rows(bead)], cls="small",
                    caption="Three 96-way split-pool extensions add barcode 1, barcode 2 and barcode 3. The completed 117-nt oligo ends in the s7 capture sequence."))


def construction() -> str:
    return '<h2>Library construction</h2>\n'+protocol_graph_panels(H.protocol_graph(),cls="long")


def sequencing() -> str:
    lib = H.final_library()
    return sp.section(lib, H.SEQ_PRIMERS, heading="Sequencing",
        intro='Read 1: 50 cycles; cell-barcode Index 1: 52; sample Index 2: 10; Read 2: 50. All four use standard Nextera sequencing primers.',
        read_lengths={"Read 1":50,"Index 1 (i7)":52,"Index 2 (i5)":10,"Read 2":50},
        show_read_content=True)


def render() -> str:
    return "\n".join([
        head("HyDrop-ATAC library chemistry"), '<div class="wrap" data-protocol-complete="true">',
        '<h1>HyDrop-ATAC &mdash; droplet single-cell chromatin accessibility</h1>',
        info('Defining source: <a href="https://doi.org/10.7554/eLife.73971">De Rop et al., <i>eLife</i> (2022)</a>.'),
        oligos(), construction(), sequencing(), '</div>'
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
