"""Construct model for single-cell reduced-representation bisulfite sequencing.

The wet-lab order is from Guo et al. (Genome Research 2013). The paper identifies
standard premethylated indexed Illumina adapters; their single-index sequences are
published in Illumina's authoritative adapter-sequence document.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from rrbs import (CELL_INDEX, INDEXED_FIXED_5, INDEX_NT, MSPI_OVERHANG,
                  MSPI_SITE, SEQ_PRIMERS, UNIVERSAL, a_tailed_fragment_scene,
                  adapter_scene, digested_fragment_scene, final_library,
                  indexed_segments, universal_segments)
from adapter_ligation import reported_adapter_end
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1101/gr.161679.113"
    graph = ProtocolGraph("scRRBS", authority=source)
    genome = graph.input("genome", "single-cell genomic DNA",
                         MaterialKind.BIOLOGICAL, primary=True)
    mspi = graph.input("mspi", "MspI restriction enzyme", MaterialKind.CATALYST)
    digested = graph.product("mspi_fragments", "MspI-digested genomic fragments",
                             MaterialKind.MOLECULAR, state=digested_fragment_scene())
    graph.react("mspi_digest", "Digest genomic DNA at C^CGG",
                StepKind.DIGESTION,
                ((genome, InputRole.TEMPLATE), (mspi, InputRole.CATALYST)), (digested,))
    klenow = graph.input("klenow", "Klenow exo- end-repair/dA-tailing enzyme",
                         MaterialKind.CATALYST)
    prepared = graph.product("a_tailed_fragments", "filled and dA-tailed MspI fragments",
                             MaterialKind.MOLECULAR, state=a_tailed_fragment_scene())
    graph.react("end_prepare", "Fill recessed ends and add one 3'-A",
                StepKind.END_REPAIR,
                ((digested, InputRole.TEMPLATE), (klenow, InputRole.CATALYST)),
                (prepared,))
    adapters = graph.input("methylated_adapters", "premethylated indexed TruSeq adapters",
                           MaterialKind.ADAPTER, state=adapter_scene())
    ligase = graph.input("ligase", "DNA ligase", MaterialKind.CATALYST)
    ligated = graph.product("ligated_fragments", "adapter-ligated genomic fragments",
                            MaterialKind.MOLECULAR)
    graph.react("adapter_ligation", "Ligate premethylated indexed Y-adapters",
                StepKind.LIGATION,
                ((prepared, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,),
                geometry=reported_adapter_end(
                    "premethylated TruSeq Y adapter", kind="cohesive",
                    polarity="3-prime", terminal_sequence="T",
                    substrate_scenes=(a_tailed_fragment_scene(),),
                    statement="The standard premethylated TruSeq adapter presents one 3'-dT to each dA-tailed insert end."))
    bisulfite = graph.input("bisulfite", "bisulfite-conversion conditions",
                            MaterialKind.CONDITION)
    converted = graph.product("converted_library", "bisulfite-converted ligation products",
                              MaterialKind.MOLECULAR)
    graph.react("bisulfite_conversion", "Bisulfite-convert adapter-ligated molecules",
                StepKind.CHEMICAL_CONVERSION,
                ((ligated, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)),
                (converted,))
    qp1 = graph.input("qp1", "QP1 / P5 PCR primer", MaterialKind.OLIGO)
    qp2 = graph.input("qp2", "QP2 / P7 PCR primer", MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "uracil-tolerant PCR polymerase",
                             MaterialKind.CATALYST)
    final = final_library()
    library = graph.product("library", final.name, MaterialKind.LIBRARY, state=final)
    graph.react("library_pcr", "Amplify converted molecules and complete flow-cell arms",
                StepKind.PCR,
                ((converted, InputRole.TEMPLATE), (qp1, InputRole.STRUCTURAL),
                 (qp2, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="The converted pre-PCR intermediate is not represented base-for-base; the final canonical arms are exact."))
    graph.sequence(library, SequencingPlan(
        library=final, primers=tuple(SEQ_PRIMERS),
        acquired_roles=("Read 1", "Index 1 (i7)", "Read 2"),
        read_lengths={"Index 1 (i7)": INDEX_NT},
        read_uses={"Read 1": "include", "Index 1 (i7)": "demultiplex",
                   "Read 2": "include"}))
    graph.declare_reviewed_scope((
        "genome", "mspi", "mspi_digest", "klenow", "end_prepare",
        "methylated_adapters", "ligase", "adapter_ligation", "bisulfite",
        "bisulfite_conversion", "qp1", "qp2", "pcr_polymerase", "library_pcr",
        "library"), source)
    return graph.require_publishable()
