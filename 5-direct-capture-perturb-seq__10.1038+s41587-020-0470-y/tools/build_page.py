#!/usr/bin/env python3
"""Build the 5-prime direct-capture Perturb-seq schematic."""
from __future__ import annotations
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import direct_capture as C
from chemdraw import Scene, annotation_rows, oligo, panel
from page import head, info
from protocol_views import assemble_page, WorkflowSection, PrimerBindingSection

OUT = HERE.parent / "5-direct-capture-perturb-seq.html"


def render() -> str:
    cdna = C.amplified_guide_cdna()
    graph = C.protocol_graph()
    return assemble_page(graph, [
        head("5′ direct-capture Perturb-seq"), '<div class="wrap" data-protocol-complete="true">',
        "<h1>5′ direct-capture Perturb-seq</h1>",
        '<p class="research-notes"><a href="01_5p-direct-capture-perturb-seq.html">Research notes</a></p>',
        info('Defining source: <a href="https://doi.org/10.1038/s41587-020-0470-y">Replogle et al., <i>Nature Biotechnology</i> (2020)</a>. This is the standard, unmodified sgRNA-CR1 branch of the 5′ protocol.'),
        "<h2>Key oligos</h2>",
        f'<seq>{oligo("oJR160", C.ojr160_segments())}</seq>',
        f'<seq>{oligo("oJR163 (forward)", C.ojr163_segments())}\n\n{oligo("oJR165 (indexed reverse)", C.ojr165_segments())}</seq>',
        WorkflowSection("Reaction workflow"),
        "<h2>Final guide library</h2>",
        panel(C.final_scene().rows(), cls="small",
              caption=f"Final guide library: {len(C.final_library())} bp for the printed AGGAGTCC example index (reported as approximately 250 bp)."),
        PrimerBindingSection("library", "Sequencing", "Sequencing",
                             intro_html=("Read 1 records cell barcode and UMI. Read 2 "
                                         "crosses the constant region before reaching the "
                                         "protospacer; the protocol recommends 98 Read 2 cycles.")),
        "</div>",
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
