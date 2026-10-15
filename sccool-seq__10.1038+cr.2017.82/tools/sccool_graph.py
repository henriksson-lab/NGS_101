"""Source-owned reaction graph for scCOOL-seq."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "scnome-seq__10.7554+eLife.23203" / "tools"))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import nomecool as N
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1038/cr.2017.82"
    graph = ProtocolGraph(N.SCCOOL, authority=source)
    converted = N.mark_and_convert_graph(graph)

    primer1 = graph.input("random_primer1", "biotinylated random primer 1",
                          MaterialKind.OLIGO, state=tuple(N.primer1_segments()))
    polymerase = graph.input("extension_polymerase", "random-priming polymerase",
                             MaterialKind.CATALYST)
    first = graph.product("first_strands", "biotinylated first-strand products",
                          MaterialKind.MOLECULAR)
    graph.react("first_priming", "Random-prime the first strand", StepKind.EXTENSION,
                ((converted, InputRole.TEMPLATE), (primer1, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (first,))

    beads = graph.input("streptavidin_beads", "streptavidin capture beads",
                        MaterialKind.SELECTION)
    captured = graph.product("captured_strands", "bead-bound first strands",
                             MaterialKind.MOLECULAR)
    waste = graph.product("uncaptured_material", "uncaptured material", MaterialKind.WASTE)
    graph.react("capture", "Capture biotinylated first strands and wash",
                StepKind.AFFINITY_CAPTURE,
                ((first, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (captured, waste))
    graph.discard(waste, "unbound material and removed converted templates are washed away")

    primer2 = graph.input("random_primer2", "random primer 2", MaterialKind.OLIGO,
                          state=tuple(N.primer2_segments()))
    template = graph.product("dual_handle_products", "dual-handle second-strand products",
                             MaterialKind.MOLECULAR)
    graph.react("second_priming", "Random-prime the captured second strand",
                StepKind.EXTENSION,
                ((captured, InputRole.TEMPLATE), (primer2, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (template,))

    forward = graph.input("forward_pe1", "Forward PE1.0 PCR primer", MaterialKind.OLIGO,
                          inferred=True, source=source + " (sequence inferred from reagent name)")
    reverse = graph.input("indexed_reverse", "NEB indexed reverse PCR primer",
                          MaterialKind.OLIGO, inferred=True,
                          source=source + " (exact primer sequence unavailable)")
    reviewed = ("genome", "mcvi", "gpc_marking", "bisulfite", "bisulfite_conversion",
                "random_primer1", "extension_polymerase", "first_priming",
                "streptavidin_beads", "capture", "random_primer2", "second_priming",
                "forward_pe1", "indexed_reverse")
    return N.finish_indexed_graph(graph, template, forward, reverse,
                                  source=source, reviewed=reviewed)
