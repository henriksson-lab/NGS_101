#!/usr/bin/env python3
"""Build the FLASH-seq protocol schematic."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHARED = HERE.parents[1] / "smart-seq__10.1038+nbt.2282" / "tools"
sys.path.insert(0, str(SHARED))

import page_parts as P  # noqa: E402
from protocol_graph import ProtocolGraph  # noqa: E402

OUT = HERE.parent / "flash-seq.html"


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: doi:10.1038/s41587-022-01312-3"
    graph = ProtocolGraph("FLASH-seq", authority=source)
    mrna = P.input_mrna(graph)
    first, enzyme, switch_scene, template = P.reverse_transcribe(
        graph, "FLASH-seq", mrna,
        rt_label="biotinylated anchored oligo-dTVN with SMART handle")
    switched = P.template_switch(
        graph, first, enzyme, switch_scene, template,
        tso_label="FLASH-seq template-switch oligo",
        step_label="Template-switch at the transcript 5-prime end")
    amplified = P.preamplify(
        graph, switched, template, family="ss2",
        statement="The source-defined SMART handle occurs in opposite orientations at the two cDNA ends, so one ISPCR reagent supplies both inward sites.")
    tagged = P.tagment(graph, amplified)
    filled = P.fill_tn5_gaps(graph, tagged)
    P.index_pcr(graph, filled, P.s2.final_library(), tuple(P.s2.SEQ_PRIMERS),
                name="FLASH-seq")
    return P.review_smart_graph(graph, source, family="ss2")


def main() -> None:
    OUT.write_text(P.render_page("FLASH-seq", graph_factory=protocol_graph), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
