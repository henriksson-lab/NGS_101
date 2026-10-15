"""MARS-seq2.0 reaction graph, owned by its defining protocol."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "mars-seq__10.1126+science.1247651" / "tools"), str(HERE.parents[1] / "lib")]
import marsseq as M
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind)
from rna_amplification import barcoded_ivt_front


def protocol_graph() -> ProtocolGraph:
    protocol = "MARS-seq2.0"
    source = "https://doi.org/10.1038/s41596-019-0164-4"
    graph = ProtocolGraph(protocol, authority=source)
    primary = graph.input("mrna", "polyadenylated RNA from one well",
                          MaterialKind.BIOLOGICAL, state=M.rt_scene(protocol),
                          primary=True)
    arna, reviewed = barcoded_ivt_front(
        graph, primary, rt_primer_label="MARS-seq2.0 RT1",
        rt_state=M.rt_scene(protocol),
        first_state=M.rt_scene(protocol),
        duplex_state=M.ds_cdna(protocol), arna_state=M.amplified_rna(protocol))
    reviewed.insert(0, "mrna")

    fragmentation = graph.input("fragmentation_condition",
                                "RNA fragmentation condition",
                                MaterialKind.CONDITION)
    fragmented = graph.product("fragmented_arna", "fragmented amplified RNA",
                               MaterialKind.MOLECULAR,
                               state=M.amplified_rna(protocol))
    graph.react("fragment_arna", "Fragment amplified RNA",
                StepKind.FRAGMENTATION,
                ((arna, InputRole.TEMPLATE),
                 (fragmentation, InputRole.CONDITION)), (fragmented,))

    adapter = graph.input("pool_adapter", "pool-barcoded RNA ligation adapter",
                          MaterialKind.ADAPTER, state=M.mars2_ligation_adapter())
    ligase = graph.input("rna_ligase", "RNA ligase", MaterialKind.CATALYST)
    ligated = graph.product("ligated_arna", "pool-adapter-ligated amplified RNA",
                            MaterialKind.MOLECULAR,
                            state=M.ligated_arna(protocol))
    geometry = AdapterLigationGeometry.unavailable(
        "The source prints the adapter sequence and terminal modifications, but the "
        "compact model does not establish both ligation-facing RNA end geometries.")
    graph.react("ligate_pool_adapter", "Ligate the pool-barcoded RNA adapter",
                StepKind.LIGATION,
                ((fragmented, InputRole.TEMPLATE),
                 (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,), geometry=geometry)

    rt2 = graph.input(
        "rt2", "RT2 primer (sequence not disclosed)", MaterialKind.OLIGO,
        inferred=True, source=source + " (RT2 sequence not disclosed)")
    reverse_transcriptase = graph.input("rt2_reverse_transcriptase",
                                       "RT2 reverse transcriptase",
                                       MaterialKind.CATALYST)
    final = M.final_library(protocol)
    pcr_template = graph.product(
        "pcr_template", "second-strand library-PCR template",
        MaterialKind.MOLECULAR, state=final, inferred=True,
        source=source + " (RT2 and PCR-primer sequences not disclosed)")
    graph.react("second_reverse_transcription",
                "Reverse-transcribe the adapter-ligated amplified RNA with RT2",
                StepKind.REVERSE_TRANSCRIPTION,
                ((ligated, InputRole.TEMPLATE), (rt2, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (pcr_template,))

    primer_source = source + " (library-PCR primer sequences not disclosed)"
    forward = graph.input("pcr_forward", "library-PCR forward primer",
                          MaterialKind.OLIGO, inferred=True, source=primer_source)
    reverse = graph.input("pcr_reverse", "library-PCR reverse primer",
                          MaterialKind.OLIGO, inferred=True, source=primer_source)
    pcr_geometry = PCRGeometry.unavailable(
        final, reason="The defining source does not disclose both PCR primer sequences.")
    polymerase = graph.input("pcr_polymerase", "PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(
        "library", f"{protocol} sequencing library", MaterialKind.LIBRARY,
        state=final, inferred=True,
        source=source + " (outer PCR arms not disclosed)")
    graph.react("library_pcr", "Amplify the final sequencing library",
                StepKind.PCR,
                ((pcr_template, InputRole.TEMPLATE),
                 (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library,),
                geometry=pcr_geometry)
    graph.sequence(library, SequencingPlan(
        limitation=("The source declares Read 1 and Read 2 contents, but does not "
                    "disclose exact sequencing-primer sites.")))
    reviewed.extend((
        "fragmentation_condition", "fragment_arna", "pool_adapter", "rna_ligase",
        "ligate_pool_adapter", "rt2", "rt2_reverse_transcriptase",
        "second_reverse_transcription", "pcr_forward", "pcr_reverse",
        "pcr_polymerase", "library_pcr", "library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
