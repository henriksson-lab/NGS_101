#!/usr/bin/env python3
"""Build the diagram-centric page for transcriptome-wide RNA FISSEQ."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "lib"))

import fisseq as F
import seqprimers as sp
from chemdraw import panel, protocol_graph_panels
from page import head, info, table

OUT = HERE.parent / "fisseq.html"


def sequencing() -> str:
    rows = [(name, f"5′-p {seq}-3′", len(seq)) for name, seq in F.SEQUENCING_PRIMERS]
    return f'''<h2>Sequencing primers on the finished rolony</h2>
{panel(F.sequencing_scene().rows(), cls="long", caption="All five sliding primers are shown on the middle repeat. The three primer-N sites in this three-repeat drawing are found computationally; every physical repeat carries the same nested sites.")}
{info("SOLiD sequencing-by-ligation uses primer N, then four primers shortened successively at the 5′ end. The offsets make the five-base ligation stride interrogate all positions.")}
{table(("Primer", "Sequence", "Length"), rows, scroll=False)}
{sp.read_spec_table([
    sp.read_spec("Sequencing output", None,
                 "in-situ rolony cDNA interrogated by five offset SOLiD ligation primers",
                 "include")],
    "FISSEQ uses SOLiD sequencing-by-ligation rather than an R1/I1/I2/R2 run.")}
'''


def render() -> str:
    graph = F.protocol_graph()
    return "\n".join([head("FISSEQ library chemistry"),
        '<div class="wrap" data-protocol-complete="true">',
        '<h1>FISSEQ &mdash; transcriptome sequencing inside fixed cells</h1>',
        '<p class="research-notes"><a href="01_fisseq.html">Research notes</a></p>',
        info('Defining source: <a href="https://doi.org/10.1126/science.1250212">Lee et al., <i>Science</i> (2014)</a>. Exact oligos and reaction order are from the authors\' open <a href="https://doi.org/10.1038/nprot.2014.191"><i>Nature Protocols</i> procedure</a>.'),
        '<h2>Circular template and RCA primer</h2>',
        F.RCA_TOPOLOGY.folded_panel(
            "The complete cDNA is covalently closed. The protected RCA primer binds "
            "the adapter antiparallel; Phi29 extends from its 3′ end, and the resulting "
            "concatemer is shown as the next molecular state in the workflow."),
        '<h2>Reaction workflow</h2>', protocol_graph_panels(graph, cls="long"),
        sequencing(), '</div>'])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__": main()
