"""Source-owned graph for the two-level dsciATAC-seq branch."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "dscatac-seq__10.1038+s41587-019-0147-6" / "tools"))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import dscatac as D
import nextera as nx
from chemdraw import Segment
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, StepKind

ASSEMBLY = D.SOURCE + " (Methods 'Assembly of indexed Tn5 transposome complexes'; Supplementary Table 6)"
WELLS = D.SOURCE + " (Methods 'Species mixing controls': 24/48/96-plex wells, pooled and washed without Tn5)"


def protocol_graph() -> ProtocolGraph:
    graph = ProtocolGraph("dsciATAC-seq", authority=D.SOURCE)
    nuclei = graph.input("nuclei", "intact cells or nuclei, one aliquot per well", MaterialKind.BIOLOGICAL,
                         primary=True)
    r1 = graph.input("indexed_r1_adaptor",
                     "Nextera Read 1 adaptor with 6-nt barcode (one of 96, HPLC-purified)",
                     MaterialKind.ADAPTER, state=tuple(D.tn5_r1("dsciATAC-seq")), source=ASSEMBLY)
    r2 = graph.input("r2_adaptor", "Nextera Read 2 adaptor", MaterialKind.ADAPTER,
                     state=(Segment("s7", nx.S7, "s7"), Segment("ME", nx.ME, "me")), source=ASSEMBLY)
    me_comp = graph.input("me_complement",
                          "ME-complementary oligo /5Phos/CTGTCTCTTATACA/3ddC/",
                          MaterialKind.OLIGO, state=(Segment("ME complement", nx.ME_RC[:14], "me"),),
                          source=ASSEMBLY)
    anneal = graph.input("anneal_condition",
                         "1:1:2 R1:R2:ME-complement, 100 uM in 50 mM NaCl; 85 C, -1 C/min to 20 C, 2 min",
                         MaterialKind.CONDITION, source=ASSEMBLY)
    annealed = graph.product("annealed_adaptors", "annealed indexed s5 and s7 ME duplexes",
                             MaterialKind.ADAPTER, source=ASSEMBLY)
    graph.react("anneal_adaptors", "Anneal each indexed R1 adaptor and the R2 adaptor to the ME complement",
                StepKind.HYBRIDIZATION,
                ((r1, InputRole.TEMPLATE), (r2, InputRole.TEMPLATE), (me_comp, InputRole.STRUCTURAL),
                 (anneal, InputRole.CONDITION)), (annealed,), source=ASSEMBLY)
    tn5 = graph.input("tn5", "purified Tn5 (14.8 uM; adaptors diluted 1:1 in glycerol, mixed 1:1, 30 min RT)",
                      MaterialKind.CATALYST, source=ASSEMBLY)
    tagged = graph.product("tagmented_cells", "well-specific Tn5-barcoded tagmented cells",
                           MaterialKind.MOLECULAR, source=WELLS)
    graph.react("tagment", "Load Tn5 with one indexed adaptor set per well, then tagment that well's cells",
                StepKind.TAGMENTATION,
                ((nuclei, InputRole.TEMPLATE), (annealed, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged,), source=WELLS)
    wash = graph.input("wash_buffer", "tagmentation buffer without Tn5", MaterialKind.SELECTION, source=WELLS)
    pooled = graph.product("pooled_cells", "pooled, washed Tn5-barcoded cells", MaterialKind.MOLECULAR,
                           source=WELLS)
    wash_waste = graph.product("wash_waste", "wash supernatant", MaterialKind.WASTE, source=WELLS)
    graph.react("pool_wash", "Pool all wells and wash cells with tagmentation buffer without Tn5",
                StepKind.CLEANUP, ((tagged, InputRole.TEMPLATE), (wash, InputRole.SELECTION)),
                (pooled, wash_waste), source=WELLS)
    graph.discard(wash_waste, "wash supernatant is removed")
    D.droplet_library(graph, pooled, protocol="dsciATAC-seq")
    graph.declare_reviewed_scope((
        "nuclei", "indexed_r1_adaptor", "r2_adaptor", "me_complement", "anneal_condition",
        "anneal_adaptors", "tn5", "tagment", "wash_buffer", "pool_wash", *D.DROPLET_SCOPE), D.SOURCE)
    return graph.require_publishable()
