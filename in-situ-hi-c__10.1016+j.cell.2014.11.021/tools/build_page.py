#!/usr/bin/env python3
"""Build the diagram-centric page for Rao et al. in situ Hi-C."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "lib"))

import insitu_hic as H
import seqprimers as sp
from chemdraw import Scene, annotation_rows, junction_row, panel, strand_row
from page import head, info

OUT = HERE.parent / "in-situ-hi-c.html"


def digest() -> str:
    return f'''<h2>Restriction digest inside intact nuclei</h2>
<h3>(1) MboI cuts the two strands on opposite sides of GATC</h3>
{panel(H.cleavage_scene().rows(), cls="small", caption="The cut coordinates belong to MboI: ^GATC on the top strand and CTAG^ on the bottom strand.")}
<h3>(2) The two products carry complementary 5′ GATC overhangs</h3>
{panel(H.cut_left_scene().rows(), cls="small", caption="Left-hand restriction fragment.")}
{panel(H.cut_right_scene().rows(), cls="small", caption="Right-hand restriction fragment.")}
'''


def fill_and_ligate() -> str:
    junction = H.contact_junction()
    rows = [*Scene.duplex(list(junction), label="contact product").rows(),
            junction_row(junction, "left filled suffix", "right filled prefix",
                         "proximity ligation"), *annotation_rows(junction)]
    return f'''<h2>Mark the ends, then proximity-ligate</h2>
<h3>(3) Klenow fills each recessed 3′ end</h3>
{panel(H.filled_left_scene().rows(), cls="small", caption="The left fragment receives GATC on its top strand; its A is biotin-14-dATP.")}
{panel(H.filled_right_scene().rows(), cls="small", caption="The right fragment receives CTAG on its bottom strand; its A is biotin-14-dATP.")}
<h3>(4) Blunt ends that are close in the nucleus ligate</h3>
{panel(rows, cls="small", caption="A contact between two filled MboI ends creates the derived GATCGATC junction. The two biotins lie on opposite strands; ** marks the exact ligation boundary.")}
'''


def library() -> str:
    insert, lib = H.sheared_insert(), H.final_library()
    return f'''<h2>Enrich contact junctions and make the sequencing library</h2>
<h3>(5) Reverse crosslinks, shear, and capture biotin</h3>
{panel([strand_row(insert), *annotation_rows(insert)], cls="small", caption="DNA is sheared to 300–500 bp. Streptavidin retains fragments containing the internal biotin-marked contact junction.")}
<h3>(6) Repair ends, add dA, ligate an indexed Illumina adapter, and PCR</h3>
{panel([*Scene.duplex(list(lib), label="library").rows(), *annotation_rows(lib)], cls="small", caption="PCR-completed paired-end library. Dotted adapter regions are the canonical TruSeq structure: the Rao protocol names an Illumina indexed adapter but does not print its bases.")}
'''


def sequencing() -> str:
    return sp.section(H.final_library(), H.SEQ_PRIMERS, heading="Sequencing",
        intro="Paired reads begin in the two genomic loci at opposite ends of the selected molecule. The proximity-ligation junction normally lies inside the insert; the two reads are mapped separately to recover the contacting loci.",
        required_roles=("Read 1", "Index 1 (i7)", "Read 2"))


def render() -> str:
    return "\n".join([head("in situ Hi-C library chemistry"), '<div class="wrap">',
        '<h1>in situ Hi-C &mdash; MboI proximity ligation in intact nuclei</h1>',
        '<p class="research-notes"><a href="01_in-situ-hi-c.html">Research notes</a></p>',
        info('Defining source: <a href="https://doi.org/10.1016/j.cell.2014.11.021">Rao et al., <i>Cell</i> (2014)</a>.'),
        info('Crosslink chromatin, digest with MboI, fill the cohesive ends with biotin-dATP, and ligate while the nuclei remain intact.'),
        digest(), fill_and_ligate(), library(), sequencing(), '</div>'])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
