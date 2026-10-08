#!/usr/bin/env python3
"""Build the PacBio SMRTbell prep kit 3.0 schematic."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "lib"))

import smrtbell as S
from chemdraw import Row, annotation_rows, panel, strand_row
from dumbbell import dumbbell_rows
from page import head, info

OUT = HERE.parent / "smrtbell-prep-kit-3-0.html"


def render() -> str:
    cycle = S.polymerase_cycle()
    primer_rows = [
        Row(chunks=[("sequencing primer + polymerase  --->  [left hairpin primer site]", "r2", False)]),
        *dumbbell_rows(S.SMRTBELL),
        Row(chunks=[("                                      [right hairpin primer site]  <---  sequencing primer + polymerase", "r2", False)]),
    ]
    return "\n".join([head("PacBio SMRTbell prep kit 3.0"), '<div class="wrap">',
        '<h1>PacBio SMRTbell prep kit 3.0 &mdash; closed dumbbell libraries</h1>',
        '<p class="research-notes"><a href="01_smrtbell-prep-kit-3-0.html">Research notes</a></p>',
        info('Commercial protocol: PacBio WGS/metagenome procedure <a href="https://www.pacb.com/wp-content/uploads/Procedure-checklist-Preparing-whole-genome-and-metagenome-libraries-using-SMRTbell-prep-kit-3.0.pdf">102-166-600 Rev09</a>.'),
        '<h2>Prepare the insert and close both ends</h2>',
        '<h3>(1) Shear, repair damage and ends, phosphorylate, and dA-tail</h3>',
        panel(S.end_prep_scene().rows(), cls="small", caption="The recommended WGS insert is 15–20 kb. Repair produces 5′ phosphates and one 3′ dA on both ends."),
        '<h3>(2) Ligate a T-overhang hairpin adapter to each end</h3>',
        panel(dumbbell_rows(S.SMRTBELL), cls="small", caption="The ds insert plus two hairpins is one covalently closed ssDNA path. ** marks the four sealed insert–adapter junctions. Adapter bases are proprietary, so the loops are structural labels rather than invented sequence."),
        '<h3>(3) Nuclease removes every molecule that still has a free end</h3>',
        info('The 37 °C nuclease treatment removes unligated DNA and leftover adapters. A completely closed SMRTbell has no free 5′ or 3′ end and survives.'),
        '<h2>Sequencing-primer and polymerase binding</h2>',
        panel(primer_rows, cls="long", caption="Both identical hairpin adapters contain a standard sequencing-primer site. The current adapter and primer sequences are not published; both possible sites are shown without fabricated bases."),
        '<h3>One productive polymerase trajectory</h3>',
        panel([strand_row(cycle), *annotation_rows(cycle)], cls="long", caption="The polymerase reads insert strand 1, traverses the opposite hairpin, reads the reverse-complement strand, traverses the first hairpin, and continues for repeated passes. This order is generated from the dumbbell topology."),
        '</div>'])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
