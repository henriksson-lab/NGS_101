#!/usr/bin/env python3
"""Build the NEBNext Ultra II USER-adaptor schematic."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "lib"))

import nebnext_user as N
import seqprimers as sp
from chemdraw import Row, Scene, annotation_rows, oligo, panel, strand_row
from dumbbell import dumbbell_rows
from page import head, info

OUT = HERE.parent / "nebnext-ultra-ii-user.html"


def render() -> str:
    hp, lib = N.hairpin_oligo(), N.final_library()
    opened = [
        Row(chunks=[("5'-p ", None, False), (N.OPENED.left_arm, "r2", False),
                    (" -3'-P   [one-nucleotide gap]   P-5'- ", None, False),
                    (N.OPENED.right_arm, "r1", False), (" -3'", None, False)])
    ]
    return "\n".join([head("NEBNext Ultra II DNA library prep — USER adaptor"), '<div class="wrap">',
        '<h1>NEBNext Ultra II DNA Library Prep &mdash; USER-cleavable adaptor</h1>',
        '<p class="research-notes"><a href="01_nebnext-ultra-ii-user.html">Research notes</a></p>',
        info('Commercial workflow: NEBNext Ultra II DNA Library Prep E7645/E7103 with the public E7600 hairpin adaptor and dual-index primers.'),
        '<h2>End prep and hairpin ligation</h2>',
        '<h3>(1) Repair ends, phosphorylate, and add one 3′ dA</h3>',
        panel(N.end_prep_scene().rows(), cls="small", caption="End prep creates the complementary dA overhangs and ligatable 5′ phosphates."),
        '<h3>(2) Ligate the single-oligo hairpin adaptor to both ends</h3>',
        '<seq>' + oligo("NEBNext Adaptor", list(hp), mods="/5Phos/", three="-3′; * = phosphorothioate bond before the terminal T") + '</seq>',
        panel([strand_row(hp), *annotation_rows(hp)], cls="long", caption="One 65-nt oligo: a verified 12-bp stem, 40-nt loop with dU, and phosphorothioate-protected 3′ dT overhang."),
        panel(dumbbell_rows(N.PRE_USER), cls="small", caption="Before USER, two ligated hairpins make a covalently closed dumbbell. ** marks the sealed insert–adaptor junctions."),
        '<h2>USER opens both loops</h2>',
        '<h3>(3) UDG removes dU; Endonuclease VIII cleaves the abasic backbone</h3>',
        panel(opened, cls="long", caption="One opened adaptor is shown. The same computed split occurs at both ends, leaving a one-nucleotide gap with 3′-phosphate and 5′-phosphate termini."),
        '<h2>Dual-index PCR and sequencing</h2>',
        '<h3>(4) i501 and i701 complete the truncated adaptor</h3>',
        info('At least three PCR cycles are required. This concrete E7600 example adds i5 TATAGCCT and i7 CGAGTAAT; the final i7 index read is its reverse complement, ATTACTCG.'),
        panel([*Scene.duplex(list(lib), label="library").rows(), *annotation_rows(lib)], cls="long", caption="PCR-completed four-read Illumina library."),
        sp.section(lib, N.sequencing_primers(), intro="All four standard primers are declared by reference and placed on the completed library by sequence matching."),
        '</div>'])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
