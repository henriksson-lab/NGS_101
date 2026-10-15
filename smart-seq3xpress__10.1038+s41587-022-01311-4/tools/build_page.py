#!/usr/bin/env python3
"""Build the SMART-seq3xpress protocol schematic."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHARED = HERE.parents[1] / "smart-seq__10.1038+nbt.2282" / "tools"
sys.path.insert(0, str(SHARED))

import page_parts as P  # noqa: E402
from protocol_graph import ProtocolGraph  # noqa: E402

OUT = HERE.parent / "smart-seq3xpress.html"


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: doi:10.1038/s41587-022-01311-4"
    graph = ProtocolGraph("SMART-seq3xpress", authority=source)
    mrna = P.input_mrna(graph)
    first, enzyme, switch_scene, template = P.reverse_transcribe(
        graph, "SMART-seq3xpress", mrna, rt_label="biotinylated anchored oligo(dT) primer")
    switched = P.template_switch(
        graph, first, enzyme, switch_scene, template,
        tso_label="SMART-seq3xpress UMI template-switch oligo",
        step_label="Template-switch at the transcript 5-prime end")
    amplified = P.preamplify(
        graph, switched, template, family="ss3",
        statement="The published asymmetric forward and reverse primers face inward around the full-length template-switched cDNA.")
    tagged = P.tagment(graph, amplified)
    filled = P.fill_tn5_gaps(graph, tagged)
    P.index_pcr(graph, filled, P.ss.smartseq3_library("SMART-seq3xpress"),
                tuple(P.ss.SS3_SEQ_PRIMERS), name="SMART-seq3xpress",
                internal=P.ss.smartseq3_library("SMART-seq3xpress", five_prime=False))
    return P.review_smart_graph(graph, source, family="ss3", internal=True)


def main() -> None:
    OUT.write_text(P.render_page("SMART-seq3xpress", graph_factory=protocol_graph), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
