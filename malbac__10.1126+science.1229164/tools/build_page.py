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
from chemdraw import Scene, Segment, oligo, panel, protocol_graph_panels, strand_row
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
    return (panel([strand_row(semi, "top")], cls="small",
                  caption="A first-generation semi-amplicon carries the common handle at one end.")
            + panel(duplex.rows(), cls="small",
                    caption="Copying a semi-amplicon produces a full amplicon with complementary "
                            "27-nt ends.")
            + M.loop_topology().folded_panel(
                    "At 58 °C the complementary 27-nt ends pair intramolecularly. N marks "
                    "both the unassigned handle bases and the representative genomic tract; "
                    "the fold, lengths and end pairing are defined. The looped product is "
                    "withdrawn from further random priming."))


def render() -> str:
    oligos = "\n".join([
        oligo("MALBAC random primer", M.random_primer().segments),
        oligo("MALBAC PCR primer", M.pcr_primer().segments),
    ])
    return "\n".join([
        head("MALBAC library chemistry"),
        '<div class="wrap" data-protocol-complete="true">',
        '<h1>MALBAC &mdash; quasi-linear single-cell whole-genome amplification</h1>',
        info('Defining source: <a href="https://doi.org/10.1126/science.1229164">'
             'Zong et al., Science 2012</a>.'),
        '<div class="caveat"><b>Bounded source scope.</b> Original main-article text40 establishes a 27-nt common handle plus eight variable bases, the four temperatures and five preamplification cycles. Exact handle bases, downstream preparation and run layout remain unassigned. The identified Supplementary Online Material is not checked in this capped pass; unassigned regions are not confirmed unpublished. H denotes the handle, N the random foot and X a schematic genomic interval.</div>',
        info('<b>Published-source difference.</b> The <a href="https://pmc.ncbi.nlm.nih.gov/'
             'articles/PMC3600412/">Science author manuscript</a> describes a common '
             '27-nt handle followed by <b>8 variable bases</b>, but does not print the '
             'handle. The inventors\' <a href="https://patents.google.com/patent/'
             'US20140200146A1/en">MALBAC patent</a> prints a candidate 27-nt handle '
             'in Example II; its worked preamplification primers use <b>5 random bases '
             'followed by GGG or TTT</b>. These may describe different embodiments, so '
             'the patent sequence does not establish the Science experiment\'s primer. '
             'This page follows the Science paper\'s 27 + 8 architecture and leaves '
             'the handle bases unknown.'),
        '<h2>Primers</h2><seq>', oligos, '</seq>',
        '<h2>Step-by-step amplification</h2>',
        protocol_graph_panels(M.protocol_graph(),cls="long"),
        '<h2>Sequencing library</h2>',
        sp.unavailable_diagram(M.full_amplicon(),
            'the downstream preparation and physical sequencing-primer sites remain unassigned',
            roles=('Read 1', 'Read 2'),
            caption='Sequencing primers cannot bind the MALBAC amplification product'),
        info('MALBAC ends at amplified genomic DNA. Downstream library preparation '
             'remains a warning boundary until the original SOM and the study-specific choices are checked.'),
        sp.unresolved_read_content_table(M.full_amplicon(),
            'the defining paper leaves the downstream sequencing-library kit unspecified'),
        '</div>',
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
