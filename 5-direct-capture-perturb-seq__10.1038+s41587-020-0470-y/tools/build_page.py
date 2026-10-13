#!/usr/bin/env python3
"""Build the 5-prime direct-capture Perturb-seq schematic."""
from __future__ import annotations
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import direct_capture as C
import seqprimers as sp
from chemdraw import Scene, annotation_rows, oligo, panel, protocol_graph_panels
from page import head, info

OUT = HERE.parent / "5-direct-capture-perturb-seq.html"


def render() -> str:
    cdna = C.amplified_guide_cdna()
    graph = C.protocol_graph()
    return "\n".join([
        head("5′ direct-capture Perturb-seq"), '<div class="wrap" data-protocol-complete="true">',
        "<h1>5′ direct-capture Perturb-seq</h1>",
        '<p class="research-notes"><a href="01_5p-direct-capture-perturb-seq.html">Research notes</a></p>',
        info('Defining source: <a href="https://doi.org/10.1038/s41587-020-0470-y">Replogle et al., <i>Nature Biotechnology</i> (2020)</a>. This is the standard, unmodified sgRNA-CR1 branch of the 5′ protocol.'),
        "<h2>Key oligos</h2>",
        f'<seq>{oligo("oJR160", C.ojr160_segments())}</seq>',
        f'<seq>{oligo("oJR163 (forward)", C.ojr163_segments())}\n\n{oligo("oJR165 (indexed reverse)", C.ojr165_segments())}</seq>',
        "<h2>Reaction workflow</h2>", protocol_graph_panels(graph, cls="long"),
        "<h2>Final guide library</h2>",
        panel(C.final_scene().rows(), cls="small",
              caption=f"Final guide library: {len(C.final_library())} bp for the printed AGGAGTCC example index (reported as approximately 250 bp)."),
        sp.section(C.final_library(), C.sequencing_primers(), heading="Sequencing",
                   intro="Read 1 records cell barcode and UMI. Read 2 crosses the constant region before reaching the protospacer; the protocol recommends 98 Read 2 cycles.",
                   required_roles=("Read 1", "Index 1 (i7)", "Read 2"),
                   read_lengths={"Read 1": 26, "Index 1 (i7)": 8, "Read 2": 98},
                   show_read_content=True),
        "</div>",
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
