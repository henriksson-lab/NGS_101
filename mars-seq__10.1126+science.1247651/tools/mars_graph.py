"""MARS-seq source-owned reaction graphs for the 2014 and 2.0 methods."""
from __future__ import annotations

import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]

import marsseq as M
from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Segment
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from rna_amplification import barcoded_ivt_front


def seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


PAPERS = {
    "MARS-seq": ("10.1126/science.1247651", "Science", "2014"),
    "MARS-seq2.0": ("10.1038/s41596-019-0164-4", "Nature Protocols", "2019"),
}

# The primary 2014 Supplementary Table S7 prints each full primer in a PDF table.
# Text extraction splits the 5′ and 3′ sequence fragments onto separate lines with
# “PCR forward”/“PCR reverse” between them, so whole-string source matching cannot
# see the ordered oligo. The original table was read directly for these two inputs.
SOURCE_OLIGO_EXCEPTIONS = {
    "pcr_forward": "Jaitin 2014 Supplementary Table S7, P5_Rd1: PDF extraction inserts the role label between its two sequence lines.",
    "pcr_reverse": "Jaitin 2014 Supplementary Table S7, P7_Rd2: PDF extraction inserts the role label between its two sequence lines.",
}


def graph_for(protocol: str) -> ProtocolGraph:
    if protocol != "MARS-seq":
        raise ValueError(f"unknown MARS-seq protocol: {protocol!r}")
    source = f"https://doi.org/{PAPERS[protocol][0]}"
    graph = ProtocolGraph(protocol, authority=source)
    primary = graph.input("mrna", "polyadenylated RNA from one well",
                          MaterialKind.BIOLOGICAL, state=M.rt_scene(protocol),
                          primary=True)
    arna, reviewed = barcoded_ivt_front(
        graph, primary, rt_primer_label=f"{protocol} RT1",
        rt_state=M.MARS1_RT1,
        first_state=M.rt_scene(protocol),
        duplex_state=M.ds_cdna(protocol), arna_state=M.amplified_rna(protocol),
        pool_before_second=True, exonuclease_when="after_pool",
        final_pool_after_exonuclease=True)
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

    adapter_state = M.POOL_ORDERED + "N" * M.MARS1_POOL_RANDOM_NT + M.LIG_CONSTANT
    adapter = graph.input("pool_adapter", "pool-barcoded RNA ligation adapter",
                          MaterialKind.ADAPTER, state=adapter_state)
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

    rt2 = graph.input("rt2", "RT2 primer", MaterialKind.OLIGO,
                      state=Construct([seg("RT2", M.MARS1_RT2_PRINTED)]))
    reverse_transcriptase = graph.input("rt2_reverse_transcriptase",
                                       "RT2 reverse transcriptase",
                                       MaterialKind.CATALYST)
    final = M.final_library(protocol)
    pcr_template = graph.product(
        "pcr_template", "second-strand library-PCR template",
        MaterialKind.MOLECULAR, state=final, source=source, inferred=True)
    graph.react("second_reverse_transcription",
                "INFERRED — reverse-transcribe adapter-ligated RNA with printed RT2; exact landing unresolved",
                StepKind.REVERSE_TRANSCRIPTION,
                ((ligated, InputRole.TEMPLATE), (rt2, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (pcr_template,))

    f_primer = PCRPrimer("P5_Rd1 PCR", M.MARS1_P5_PCR)
    r_primer = PCRPrimer("P7_Rd2 PCR", M.MARS1_P7_PCR)
    forward = graph.input("pcr_forward", f_primer.name, MaterialKind.OLIGO,
                          state=Construct([seg(f_primer.name, f_primer.sequence)]))
    reverse = graph.input("pcr_reverse", r_primer.name, MaterialKind.OLIGO,
                          state=Construct([seg(r_primer.name, r_primer.sequence)]))
    pcr_geometry = PCRGeometry.exact(
        final, f_primer, r_primer,
        statement="Supplementary Tables S7–S8 print both library-PCR primers.")
    polymerase = graph.input("pcr_polymerase", "PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(
        "library", f"{protocol} sequencing library", MaterialKind.LIBRARY,
        state=final, source=source, inferred=True)
    graph.react("library_pcr", "Amplify the final sequencing library",
                StepKind.PCR,
                ((pcr_template, InputRole.TEMPLATE),
                 (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library,),
                geometry=pcr_geometry)
    graph.sequence(library, SequencingPlan(
        library=final, primers=M.SEQ_PRIMERS,
        acquired_roles=("Read 1", "Read 2")))
    reviewed.extend((
        "fragmentation_condition", "fragment_arna", "pool_adapter", "rna_ligase",
        "ligate_pool_adapter", "rt2", "rt2_reverse_transcriptase",
        "second_reverse_transcription", "pcr_forward", "pcr_reverse",
        "pcr_polymerase", "library_pcr", "library"))
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()



def protocol_graph() -> ProtocolGraph:
    """The named 2014 protocol, for source and chemistry audits."""
    return graph_for("MARS-seq")
