#!/usr/bin/env python3
"""Build the original Micro-C protocol schematic."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "lib"))

import micro_c as M
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from chemdraw import (Scene, annotation_rows, junction_row, panel,
                      protocol_graph_panels, strand_row)
from page import head, info, table

OUT = HERE.parent / "micro-c.html"


def render():
    junction = M.contact_junction(); insert = M.selected_insert(); lib = M.final_library()
    graph = M.protocol_graph()
    junction_rows = [*Scene.duplex(list(junction), label="contact product").rows(),
                     junction_row(junction, "repaired end A", "repaired end B",
                                  "proximity ligation"), *annotation_rows(junction)]
    return "\n".join([
        head("Micro-C library chemistry"), '<div class="wrap" data-protocol-complete="true">',
        '<h1>Micro-C &mdash; nucleosome-resolution proximity ligation</h1>',
        '<p class="research-notes"><a href="01_micro-c.html">Research notes</a></p>',
        info('Defining source: <a href="https://doi.org/10.1016/j.cell.2015.05.048">Hsieh et al., <i>Cell</i> (2015)</a>.'),
        info('Bounded main-source check: Antarctic phosphatase precedes T4 polymerase processing. Digestion products have unresolved recovery; digestion is not proof of complete physical removal. Supporting supplements and full molecular assembly remain deferred.'),
        '<h2>Reaction workflow</h2>', protocol_graph_panels(graph, cls="long"),
        '<h2>Final sequencing library</h2>',
        panel([*Scene.duplex(list(lib), label="library").rows(), *annotation_rows(lib)],
              cls="small", caption="End repair, dA-tailing and Illumina-adapter ligation precede streptavidin capture and 12–15 PCR cycles. The paper does not identify the adapter or index architecture."),
        '<h2>Sequencing</h2><p>Paired reads begin in the two nucleosome-derived genomic ends.</p>',
        '<div class="info">Primer placement unavailable: the source specifies Illumina paired-end sequencing but not the library kit or sequencing-primer sequences.</div>',
        sp.read_spec_table([
            sp.read_spec("Read 1", M.READ_LENGTHS["Read 1"], "one nucleosome-derived genomic end", "include"),
            sp.read_spec("Read 2", M.READ_LENGTHS["Read 2"], "the other nucleosome-derived genomic end", "include")],
            "Original main Methods specifies paired50-base reads. Exact oligos and adapter assembly remain unassigned pending original supplementary source review; reaching the search cap does not establish absence."),
        '</div>'
    ])


def main():
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
