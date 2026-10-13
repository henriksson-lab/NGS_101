#!/usr/bin/env python3
"""Build the NEBNext Ultra II USER-adaptor schematic."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "lib"))

import nebnext_user as N
import seqprimers as sp
from chemdraw import (Row, Scene, annotation_rows, oligo, panel,
                      protocol_graph_panels, strand_row)
from dumbbell import dumbbell_rows
from page import head, info

OUT = HERE.parent / "nebnext-ultra-ii-user.html"


def render() -> str:
    hp, lib = N.hairpin_oligo(), N.final_library()
    graph = N.protocol_graph()
    opened = [
        Row(chunks=[("5'-p ", None, False), (N.OPENED.left_arm, "r2", False),
                    (" -3'-P   [one-nucleotide gap]   P-5'- ", None, False),
                    (N.OPENED.right_arm, "r1", False), (" -3'", None, False)])
    ]
    return "\n".join([head("NEBNext Ultra II DNA library prep — USER adaptor"), '<div class="wrap" data-protocol-complete="true">',
        '<h1>NEBNext Ultra II DNA Library Prep &mdash; USER-cleavable adaptor</h1>',
        '<p class="research-notes"><a href="01_nebnext-ultra-ii-user.html">Research notes</a></p>',
        info('Commercial workflow: NEBNext Ultra II DNA Library Prep E7645/E7103 with the public E7600 hairpin adaptor and dual-index primers.'),
        '<h2>USER-cleavable adaptor</h2>',
        '<seq>' + oligo("NEBNext Adaptor", list(hp), mods="/5Phos/", three="-3′; * = phosphorothioate bond before the terminal T") + '</seq>',
        panel([strand_row(hp), *annotation_rows(hp)], cls="long", caption="One 65-nt oligo: a verified 12-bp stem, 40-nt loop with dU, and phosphorothioate-protected 3′ dT overhang."),
        N.HAIRPIN_TOPOLOGY.folded_panel(
            "The same 65-nt adaptor folded base by base. The paired 12-nt arms form the "
            "stem; the 40-nt dU-containing interval is the loop, and the terminal 3′ T "
            "remains outside the stem as the ligating overhang."),
        '<h2>Reaction workflow</h2>', protocol_graph_panels(graph, cls="long"),
        '<h2>Dual-index library and sequencing</h2>',
        info('At least three PCR cycles are required. This concrete E7600 example adds i5 TATAGCCT and i7 CGAGTAAT; the final i7 index read is its reverse complement, ATTACTCG.'),
        panel([*Scene.duplex(list(lib), label="library").rows(), *annotation_rows(lib)], cls="long", caption="PCR-completed four-read Illumina library."),
        sp.section(lib, N.sequencing_primers(), intro="All four standard primers are declared by reference and placed on the completed library by sequence matching."),
        '</div>'])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
