#!/usr/bin/env python3
"""Build the diagram-centric scDNase-seq chemistry page."""
from __future__ import annotations

import html
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]

import scdnase_seq as S
import seqprimers as sp
from chemdraw import Scene, oligo, panel, annotation_rows, protocol_graph_panels
from page import head, info, table

OUT = HERE.parent / "scdnase-seq.html"


def sequencing() -> str:
    lib = S.final_library()
    rows = []
    for primer, hit in S.primer_landings():
        if hit is None:
            raise ValueError(f"{primer.name} has no site")
        rows.append((html.escape(primer.role), html.escape(primer.name),
                     html.escape(", ".join(hit.covers)),
                     f"<code>{html.escape(hit.reads)}</code>&hellip;"))
    return "".join([
        '<h2>Sequencing</h2>',
        info("Historical Multiplexing-kit sites are shown as a reference. The inspected Methods name HiSeq2500, but actual read and index acquisition remains unassigned."),
        sp.diagram(lib, S.SEQ_PRIMERS),
        table(("Read", "Primer", "Primer site", "First bases"), rows),
        sp.read_spec_table((sp.read_spec("Sequencing output", None, "HiSeq2500; read/index mask unassigned", "include"),)),
    ])


def render() -> str:
    lib = S.final_library()
    graph = S.protocol_graph()
    return "\n".join([
        head("scDNase-seq library chemistry"), '<div class="wrap" data-protocol-complete="true">',
        '<h1>scDNase-seq (Pico-Seq) &mdash; DNase-accessible DNA from one cell</h1>',
        info('Defining source: <a href="https://doi.org/10.1038/nature15740">Jin et al., <i>Nature</i> (2015)</a>.'),
        info('The paper’s cited library method identifies Illumina Multiplexing Sample Prep Oligonucleotide Kit ref. 1005709. The displayed obsolete six-base single-index oligos are the corresponding vendor-published set; there is no i5 index.'),
        '<h2>Key oligos</h2><seq>',
        oligo("Y-adapter, 3'-T strand", S.adapter_bottom()),
        oligo("Y-adapter, phosphorylated strand", S.adapter_top(), mods="/5Phos/"),
        oligo("PCR Primer 1.0", S.p5_primer()),
        oligo("indexed multiplexing primer", S.indexed_p7_primer()),
        '</seq>',
        info('Original Methods: DNase digestion is followed by stop buffer, proteinase K and30ng circular plasmid carrier, then phenol–chloroform extraction and ethanol precipitation with glycogen. These recovery/carrier details are not separate graph states. First PCR is6 cycles, followed by160–300bp E-gel selection,22-cycle P5/P7 PCR and a second160–300bp gel selection. The final selection is also currently summarized outside the graph.'),
        '<h2>Reaction workflow</h2>', protocol_graph_panels(graph, cls="long"),
        '<h2>Molecular details</h2>',
        '<h3>(1) Digest one lysed cell with DNase I</h3>',
        panel(S.digested_fragment_scene().rows(), cls="small",
              caption="DNase I releases short fragments from accessible genomic DNA; circular carrier DNA is added before purification."),
        '<h3>(2) End-repair and add one 3-prime A at each end</h3>',
        panel(S.a_tailed_fragment_scene().rows(), cls="small",
              caption="The dA tail provides the substrate for the kit's 3'-T adapter."),
        '<h3>(3) Ligate the forked Multiplexing adapter</h3>',
        panel(S.adapter_scene().rows(), cls="small",
              caption="The two kit oligos form a 12-bp stem and leave a 3'-T ligation overhang."),
        '<h3>(4) Index-amplify, size-select, then amplify with P5 and P7</h3>',
        panel([*Scene.duplex(list(lib), label="library").rows(), *annotation_rows(lib)],
              cls="small", caption="PCR-completed scDNase-seq library, P5 to P7'. The six-base i7 sample index identifies the cell."),
        sequencing(),
        '</div>',
    ])


def main() -> None:
    OUT.write_text(render(), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
