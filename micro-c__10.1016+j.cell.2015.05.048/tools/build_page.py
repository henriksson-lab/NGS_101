#!/usr/bin/env python3
"""Build the original Micro-C protocol schematic."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "lib"))

import micro_c as M
from chemdraw import Scene, annotation_rows, junction_row, panel, strand_row
from page import head, info, table

OUT = HERE.parent / "micro-c.html"


def render():
    junction = M.contact_junction(); insert = M.selected_insert(); lib = M.final_library()
    junction_rows = [*Scene.duplex(list(junction), label="contact product").rows(),
                     junction_row(junction, "repaired end A", "repaired end B",
                                  "proximity ligation"), *annotation_rows(junction)]
    return "\n".join([
        head("Micro-C library chemistry"), '<div class="wrap">',
        '<h1>Micro-C &mdash; nucleosome-resolution proximity ligation</h1>',
        '<p class="research-notes"><a href="01_micro-c.html">Research notes</a></p>',
        info('Defining source: <a href="https://doi.org/10.1016/j.cell.2015.05.048">Hsieh et al., <i>Cell</i> (2015)</a>.'),
        '<h2>Fragment crosslinked chromatin to mononucleosomes</h2>',
        panel(M.mnase_scene().rows(), cls="small",
              caption="MNase digests accessible linker DNA. Unlike restriction Hi-C, there is no invariant recognition-site sequence at the fragment ends."),
        '<h2>Repair and label the heterogeneous ends</h2>',
        panel(M.resected_end_scene().rows(), cls="small",
              caption="T4 DNA polymerase with ATP resects heterogeneous MNase products to 5′ single-stranded ends. One representative end is shown."),
        panel(M.repaired_end_scene().rows(), cls="small",
              caption="Adding biotin-dATP, biotin-dCTP, dGTP and dTTP fills the overhang to a labelled blunt end. Label positions are computed from the template."),
        '<h2>Proximity-ligate and remove free labelled ends</h2>',
        panel(junction_rows, cls="small",
              caption="T4 ligase joins blunt nucleosomal ends held near one another. ** marks the ligation boundary; its sequence varies with the MNase cuts."),
        table(("Molecule", "Exonuclease III result"),
              (("unligated labelled end", "terminal biotin-containing DNA is removed"),
               ("proximity-ligated product", "the now-internal labelled junction is retained")),
              scroll=False),
        '<h2>Recover contacts and build the Illumina library</h2>',
        panel([strand_row(insert), *annotation_rows(insert)], cls="small",
              caption="After crosslink reversal, the original protocol gel-selects 250–350-bp ligation products."),
        panel([*Scene.duplex(list(lib), label="library").rows(), *annotation_rows(lib)],
              cls="small", caption="End repair, dA-tailing and Illumina-adapter ligation precede streptavidin capture and 12–15 PCR cycles. The paper does not identify the adapter or index architecture."),
        '<h2>Sequencing</h2><p>Paired reads begin in the two nucleosome-derived genomic ends.</p>',
        '<div class="info">Primer placement unavailable: the source specifies Illumina paired-end sequencing but not the library kit or sequencing-primer sequences.</div>',
        '</div>'
    ])


def main():
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
