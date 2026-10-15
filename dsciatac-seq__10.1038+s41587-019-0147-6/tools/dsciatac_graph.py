"""Source-owned graph for the two-level dsciATAC-seq branch."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "dscatac-seq__10.1038+s41587-019-0147-6" / "tools"))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import dscatac as D
from protocol_graph import MaterialKind, ProtocolGraph


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1038/s41587-019-0147-6"
    graph = ProtocolGraph("dsciATAC-seq", authority=source)
    nuclei = graph.input("nuclei", "intact cells or nuclei", MaterialKind.BIOLOGICAL,
                         primary=True)
    tagged = D.tagment_cells(graph, nuclei, protocol="dsciATAC-seq")
    filled = D.fill_transposition_gaps(graph, tagged, protocol="dsciATAC-seq")
    droplets = D.partition_with_beads(graph, filled)
    D.copy_bead_barcode(graph, droplets, protocol="dsciATAC-seq", source=source)
    graph.declare_reviewed_scope((
        "nuclei", "transposome", "tn5", "tagment", "gap_fill_polymerase", "gap_fill",
        "droplet_loading", "partition", "bead_oligo", "reverse_primer", "pcr_polymerase",
        "droplet_pcr", "library"), source)
    return graph.require_publishable()
