#!/usr/bin/env python3
"""Build the diagram-centric txci-ATAC-seq page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import txci_atac as T
from chemdraw import Scene, annotation_rows, oligo, panel
from page import head, info, table

OUT = HERE.parent / "txci-atac-seq.html"


def render() -> str:
    lib = T.final_library()
    duplex = Scene.duplex(list(lib), label="library")
    landings = [(html.escape(p.role), html.escape(p.name),
                 html.escape(", ".join(hit.covers)),
                 f"<code>{html.escape(hit.reads)}</code>&hellip;")
                for p, hit in T.primer_landings(lib)]
    oligos = "\n".join([
        oligo("Tn5ME-A", T.tn5_a()),
        oligo("Tn5ME-B (one representative barcode)", T.tn5_b()),
        oligo("Tn5MErev", T.tn5_bottom(), mods="/5Phos/"),
        oligo("10x bead oligo", T.bead_oligo()),
        oligo("Short SBS", [T.seg("partial TruSeq Read 2", T.SHORT_SBS, "r2")]),
        oligo("P7.S701", T.p7_primer()),
    ])
    return "\n".join([
        head("txci-ATAC-seq library chemistry"), '<div class="wrap">',
        '<h1>txci-ATAC-seq</h1>',
        info('Defining source: <a href="https://doi.org/10.1186/s13059-023-03150-1">Zhang et al., <i>Genome Biology</i> (2024)</a>; first described in the linked preprint.'),
        '<h2>Indexed transposomes</h2>',
        panel(T.transposome_scene("A").rows(), cls="small", caption="Adaptor A supplies the s5 end recognized by the 10x bead oligo."),
        panel(T.transposome_scene("B").rows(), cls="small", caption="Adaptor B supplies one plate barcode and the partial TruSeq Read 2 site."),
        '<h2>Key oligos</h2>', '<seq>', oligos, '</seq>',
        '<h2>Library construction</h2>',
        '<h3>(1) Pre-index nuclei in a 96-well Tn5 plate</h3>',
        info('Each well has one 8-nt Tn5 barcode. Nuclei are pooled after tagmentation and overloaded into the Chromium ATAC workflow.'),
        '<h3>(2) Add GEM barcode and suppress barcode swapping</h3>',
        info('The P5&ndash;GEM barcode&ndash;s5 bead oligo primes the A end. Short SBS primes the B end, making amplification exponential within each droplet.'),
        '<h3>(3) Add the lane index and flow-cell ends</h3>',
        panel([*duplex.rows(), *annotation_rows(lib)], cls="small", caption="Final library, P5 to P7\'. Cell identity combines the GEM and Tn5 barcodes; i7 identifies the lane."),
        '<h2>Sequencing</h2>',
        info('Read 1: 51 cycles; i7: 8 usable cycles; i5: 16 cycles; Read 2: 78 cycles. Standard Illumina primer mixes are used.'),
        table(("Read", "Standard primer", "Primer site", "First bases"), landings),
        '<h3>Read 2 before genomic DNA</h3>',
        table(("Cycles", "Content"), T.read2_layout(), scroll=False),
        '</div>'
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
