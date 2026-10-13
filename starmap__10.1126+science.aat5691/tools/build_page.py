#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import starmap as M
import seqprimers as sp
from chemdraw import panel, protocol_graph_panels
from page import head, info
OUT=HERE.parent/"starmap.html"

def render():
    return "\n".join([head("STARmap chemistry"), '<div class="wrap" data-protocol-complete="true">',
      '<h1>STARmap</h1>', '<p class="research-notes"><a href="01_starmap.html">Research notes</a></p>',
      info('Defining source: <a href="https://doi.org/10.1126/science.aat5691">Wang et al., <i>Science</i> (2018)</a>. Symbolic bases mark panel-specific sequences not exposed by the available defining-paper tables.'),
      '<h2>Paired SNAIL recognition</h2>', panel(M.snail_scene().rows(), cls="long"),
      '<h2>Reaction workflow</h2>', protocol_graph_panels(M.protocol_graph(), cls="long"),
      '<h2>SEDAL interrogation</h2>', panel(M.sedal_scene().rows(), cls="long"),
      sp.read_spec_table([sp.read_spec("Sequencing output", None,
        "six microscopy cycles interrogate the 5-base target code in each anchored rolony", "include")],
        "STARmap has no Illumina R1, I1, I2 or R2 reads."), '</div>'])

if __name__=="__main__": OUT.write_text(render(), encoding="utf-8"); print(f"wrote {OUT}")
