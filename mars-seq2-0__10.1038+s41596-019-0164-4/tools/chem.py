"""MARS-seq2.0 reaction graph, owned by its defining protocol."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "mars-seq__10.1126+science.1247651" / "tools"), str(HERE.parents[1] / "lib")]
import marsseq as M
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind)
from rna_amplification import barcoded_ivt_front


# Author procedure (protocols.io 7hkhj4w, Steps 51-80), read 2026-10-17.
PROC = "https://doi.org/10.17504/protocols.io.7hkhj4w"
# Step 80: "rd1 = 75 bp, rd2 = 15 bp, no indices".  Read 2 therefore stops after the
# 7-nt well barcode and 8-nt UMI; the MOESM1 seq_batches example (3I.4P.50M /
# 7W.8R.5I) describes a different, longer configuration and is shown as such.
READ_LENGTHS = {"Read 1": 75, "Read 2": 15}


def protocol_graph() -> ProtocolGraph:
    protocol = "MARS-seq2.0"
    source = "https://doi.org/10.1038/s41596-019-0164-4"
    full_protocol = "https://doi.org/10.17504/protocols.io.7hkhj4w"
    graph = ProtocolGraph(protocol, authority=source)
    primary = graph.input("mrna", "polyadenylated RNA from one well",
                          MaterialKind.BIOLOGICAL, state=M.rt_scene(protocol),
                          primary=True)
    arna, reviewed = barcoded_ivt_front(
        graph, primary, rt_primer_label="MARS-seq2.0 RT1",
        rt_state=Construct(M.rt1_segments(protocol), name="MARS-seq2.0 RT1"),
        first_state=M.rt_scene(protocol),
        duplex_state=M.ds_cdna(protocol), arna_state=M.amplified_rna(protocol),
        pool_before_second=True, exonuclease_when="before_pool")
    reviewed.insert(0, "mrna")

    fragmentation = graph.input("fragmentation_condition",
                                "Zn2+ RNA fragmentation reagent, 3 min at 70 °C, then stop "
                                "solution (Steps 51–54; median ~300 nt)",
                                MaterialKind.CONDITION, source=PROC + " (Steps 51–55)")
    fragmented = graph.product("fragmented_arna", "fragmented amplified RNA",
                               MaterialKind.MOLECULAR,
                               state=M.amplified_rna(protocol))
    graph.react("fragment_arna", "Fragment amplified RNA",
                StepKind.FRAGMENTATION,
                ((arna, InputRole.TEMPLATE),
                 (fragmentation, InputRole.CONDITION)), (fragmented,))

    adapter = graph.input("pool_adapter", "pool-barcoded RNA ligation adapter",
                          MaterialKind.ADAPTER, state=M.mars2_ligation_adapter())
    ligase = graph.input("rna_ligase",
                         "T4 RNA ligase 1 with DMSO, ATP and PEG, 2 h at 22 °C after "
                         "adapter + RNA 3 min at 70 °C (Steps 57–60)",
                         MaterialKind.CATALYST, source=PROC + " (Steps 56–61)")
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
        "rt2", "RT2 primer", MaterialKind.OLIGO,
        state=M.MARS2_RT2, source=full_protocol + " (reagent setup, p. 29)")
    reverse_transcriptase = graph.input("rt2_reverse_transcriptase",
                                       "AffinityScript RT after RT2 annealing 3 min at "
                                       "72 °C; 2 min 42 °C, 50 min 45 °C, 5 min 85 °C "
                                       "(Steps 63–66)",
                                       MaterialKind.CATALYST,
                                       source=PROC + " (Steps 62–67)")
    final = M.final_library(protocol)
    pre_pcr = M.mars2_pcr_template()
    pcr_template = graph.product(
        "pcr_template", "second-strand library-PCR template",
        MaterialKind.MOLECULAR, state=pre_pcr,
        source=full_protocol + " (reagent setup, p. 29)")
    graph.react("second_reverse_transcription",
                "Reverse-transcribe the adapter-ligated amplified RNA with RT2",
                StepKind.REVERSE_TRANSCRIPTION,
                ((ligated, InputRole.TEMPLATE), (rt2, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (pcr_template,))

    primer_source = full_protocol + " (enrichment primers, p. 29)"
    forward = graph.input("pcr_forward", "library-PCR forward primer",
                          MaterialKind.OLIGO, state=M.MARS2_P5_PCR,
                          source=primer_source)
    reverse = graph.input("pcr_reverse", "library-PCR reverse primer",
                          MaterialKind.OLIGO, state=M.MARS2_P7_PCR,
                          source=primer_source)
    pcr_geometry = PCRGeometry.exact(
        pre_pcr, PCRPrimer("P5-Read1", M.MARS2_P5_PCR),
        PCRPrimer("P7-Read2", M.MARS2_P7_PCR),
        statement="Author protocols.io reagent setup, p. 29: RT2 and enrichment primers.")
    polymerase = graph.input("pcr_polymerase",
                             "KAPA HiFi HotStart ReadyMix; example 2 cycles at 57 °C + "
                             "12 at 65 °C annealing, typically 11–17 total (Steps 69–71)",
                             MaterialKind.CATALYST, source=PROC + " (Steps 69–72)")
    library = graph.product(
        "library", f"{protocol} sequencing library", MaterialKind.LIBRARY,
        state=final, source=full_protocol + " (enrichment primers, p. 29)")
    graph.react("library_pcr", "Amplify the final sequencing library",
                StepKind.PCR,
                ((pcr_template, InputRole.TEMPLATE),
                 (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library,),
                geometry=pcr_geometry)
    graph.sequence(library, SequencingPlan(
        acquired_roles=tuple(READ_LENGTHS), read_lengths=READ_LENGTHS,
        limitation=("The source declares read contents and the enrichment-primer "
                    "flanks; separate instrument sequencing-primer oligos are not named.")))
    reviewed.extend((
        "fragmentation_condition", "fragment_arna", "pool_adapter", "rna_ligase",
        "ligate_pool_adapter", "rt2", "rt2_reverse_transcriptase",
        "second_reverse_transcription", "pcr_forward", "pcr_reverse",
        "pcr_polymerase", "library_pcr", "library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
