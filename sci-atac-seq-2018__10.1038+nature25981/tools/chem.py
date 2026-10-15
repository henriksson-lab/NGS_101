"""Cusanovich 2018 two-round sci-ATAC source choices and graph."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] /
                       "sci-atac-seq__10.1126+science.aab1601" / "tools"))

import sciatac as S
from sciatac import indexed_transposome, pool_indexed_nuclei, seg, two_tagmented
from chemdraw import Construct
import illumina as il
import nextera as nx
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

PROTOCOL = S.SCI18
SOURCE = ("10.1038/nature25981", PROTOCOL)
PCR_INDEX_NT = 10
SOURCE_P5_1 = "AATGATACGGCGACCACCGAGATCTACACCTCCATCGAGTCGTCGGCAGCGTC"
FINAL_LIBRARY = S.two_level_library(PROTOCOL, pcr_n=PCR_INDEX_NT, inferred=False)
PRIMERS = S.TWO_LEVEL_PRIMERS
READ_ROWS = [
    ("Read 1", "50", "genomic insert"),
    ("Index 1", "8 + 27 dark + 10", "t7 barcode + PCR i7"),
    ("Index 2", "8 + 21 dark + 10", "t5 barcode + PCR i5"),
    ("Read 2", "50", "genomic insert"),
]


def protocol_graph() -> ProtocolGraph:
    protocol = PROTOCOL
    doi = SOURCE[0]
    source = "https://doi.org/" + doi
    graph = ProtocolGraph(protocol, authority=source)
    nuclei = graph.input(
        "nuclei", "accessible chromatin in nuclei", MaterialKind.BIOLOGICAL,
        state=Construct([seg("genomic DNA", "X" * 40, placeholder=True)]),
        primary=True)
    reviewed = ["nuclei"]
    transposome = graph.input(
        "indexed_transposomes", "indexed T5 and T7 transposomes",
        MaterialKind.ADAPTER, state=indexed_transposome("T5"))
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagmented = graph.product(
        "tagmented_fragments", "heterologous T5/T7 tagmentation products",
        MaterialKind.MOLECULAR, state=two_tagmented())
    tagment_waste = graph.product("tagment_waste", "same-end tagmentation products",
                                  MaterialKind.WASTE)
    graph.react("tagment", "Tagment nuclei in the first indexing plate",
                StepKind.TAGMENTATION,
                ((nuclei, InputRole.TEMPLATE),
                 (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagmented, tagment_waste))
    graph.discard(tagment_waste,
                  "same-end T5/T5 and T7/T7 products lack one terminal PCR site")
    pooled, pool_reviewed = pool_indexed_nuclei(
        graph, tagmented, prefix="round1", state=two_tagmented(),
        label="first-round indexed tagmented nuclei from other wells")
    gap_polymerase = graph.input("gap_polymerase", "gap-fill DNA polymerase",
                                 MaterialKind.CATALYST)
    final = FINAL_LIBRARY
    filled = graph.product(
        "gap_filled_fragments", "gap-filled redistributed fragments",
        MaterialKind.MOLECULAR, state=final)
    graph.react("gap_fill", "Fill the 9-nt tagmentation gaps",
                StepKind.EXTENSION,
                ((pooled, InputRole.TEMPLATE),
                 (gap_polymerase, InputRole.CATALYST)), (filled,))
    pcr_n = PCR_INDEX_NT
    p5_sequence = il.P5 + "N" * pcr_n + nx.S5
    p7_sequence = il.P7 + "N" * pcr_n + nx.S7
    p5 = graph.input(
        "p5_index_primer", "second-round P5 index primer", MaterialKind.OLIGO,
        state=Construct([seg("P5 index primer", p5_sequence, placeholder=True)]))
    p7 = graph.input(
        "p7_index_primer", "second-round P7 index primer", MaterialKind.OLIGO,
        state=Construct([seg("P7 index primer", p7_sequence, placeholder=True)]))
    polymerase = graph.input("pcr_polymerase", "PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(
        "library", f"{protocol} sequencing library", MaterialKind.LIBRARY,
        state=final)
    geometry = PCRGeometry.exact(
        final, PCRPrimer("P5 index primer", p5_sequence),
        PCRPrimer("P7 index primer", p7_sequence),
        statement="The 2018 oligo table prints both indexed PCR primers.")
    graph.react("index_pcr", "Add second-round P5/P7 indices by PCR",
                StepKind.PCR,
                ((filled, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library,),
                geometry=geometry)
    run_primers = PRIMERS
    graph.sequence(library, SequencingPlan(
        library=final, primers=run_primers,
        acquired_roles=tuple(p.role for p in run_primers)))
    reviewed.extend([
        "indexed_transposomes", "tn5", "tagment", *pool_reviewed,
        "gap_polymerase", "gap_fill", "p5_index_primer", "p7_index_primer",
        "pcr_polymerase", "index_pcr", "library",
    ])
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
