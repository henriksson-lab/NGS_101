#!/usr/bin/env python3
"""Build the MALBAC amplification schematic."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import malbac as M
import seqprimers as sp
from chemdraw import Scene, Segment, oligo, panel, strand_row
from page import head, info

OUT = HERE.parent / "malbac.html"


def random_priming() -> str:
    template = [Segment("template", "X" * 42, placeholder=True)]
    primer = list(M.random_primer().segments)
    sc = Scene()
    sc.strand("genome", template, label="genomic DNA", rev=True)
    sc.anneal("primer", primer, to="genome", pair=("N8", "template"), label="random primer")
    sc.arrow("primer", "strand-displacing extension")
    return panel(sc.rows(), cls="small",
                 caption="The N8 foot anneals at many genomic positions; the common handle "
                         "remains as a 5' overhang.")


def products() -> str:
    semi = M.semi_amplicon()
    full = M.full_amplicon()
    duplex = Scene.duplex(full.segments)
    duplex.strands["top"].label = duplex.strands["bottom"].label = ""
    folded = Scene()
    folded.strand("amplicon", full.segments, label="full amplicon")
    folded.mark("amplicon", "left handle", "complementary end 1")
    folded.mark("amplicon", "right handle'", "complementary end 2")
    return (panel([strand_row(semi, "top")], cls="small",
                  caption="A first-generation semi-amplicon carries the common handle at one end.")
            + panel(duplex.rows(), cls="small",
                    caption="Copying a semi-amplicon produces a full amplicon with complementary "
                            "27-nt ends.")
            + panel(folded.rows(), cls="small",
                    caption="At 58 °C the complementary ends pair intramolecularly. The looped "
                            "product is withdrawn from further random priming."))


def render() -> str:
    oligos = "\n".join([
        oligo("MALBAC random primer", M.random_primer().segments),
        oligo("MALBAC PCR primer", M.pcr_primer().segments),
    ])
    return "\n".join([
        head("MALBAC library chemistry"),
        '<div class="wrap">',
        '<h1>MALBAC &mdash; quasi-linear single-cell whole-genome amplification</h1>',
        info('Defining source: <a href="https://doi.org/10.1126/science.1229164">'
             'Zong et al., Science 2012</a>.'),
        '<div class="caveat"><b>Sequence availability.</b> The paper states a 27-nt common '
        'handle and an 8-nt random foot; its exact handle sequence is only in the unavailable '
        'supplement. H therefore marks a length-preserving unpublished placeholder.</div>',
        '<h2>Primers</h2><seq>', oligos, '</seq>',
        '<h2>Step-by-step amplification</h2>',
        f'<h3>(1) Lyse one cell and melt genomic DNA ({M.MELT_C}&nbsp;&deg;C)</h3>',
        f'<h3>(2) Anneal random primers ({M.ANNEAL_C}&nbsp;&deg;C) and extend '
        f'({M.EXTEND_C}&nbsp;&deg;C)</h3>', random_priming(),
        '<h3>(3) Melt to release semi-amplicons</h3>',
        f'<h3>(4) Repeat for {M.PREAMPLIFICATION_CYCLES} quasi-linear cycles; loop full '
        f'amplicons at {M.LOOP_C}&nbsp;&deg;C after each cycle</h3>', products(),
        '<h3>(5) PCR with the common 27-nt primer</h3>',
        panel(Scene.duplex(M.full_amplicon().segments).rows(), cls="small",
              caption="Amplified MALBAC product. The protocol itself adds no cell barcode, UMI "
                      "or sequencing adapter."),
        '<h2>Sequencing library</h2>',
        sp.unavailable_diagram(M.full_amplicon(),
            'a separate, unspecified library-preparation kit must first add platform adapters',
            roles=('Read 1', 'Read 2'),
            caption='Sequencing primers cannot bind the MALBAC amplification product'),
        info('MALBAC ends at amplified genomic DNA. A separate standard library-preparation kit '
             'adds the platform adapters; the accessible defining paper does not identify that kit.'),
        '</div>',
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
