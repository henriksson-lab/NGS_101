"""Source-owned microSPLiT reaction order on shared split-pool molecular states."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "split-seq__10.1126+science.aam8999" / "tools"))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import splitseq as S
from adapter_ligation import reported_adapter_end
from chemdraw import Construct, Segment
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1126/science.aba5257"
    graph = ProtocolGraph("microSPLiT split-pool RNA library", authority=source)
    first, rt, ligase = S.start_split_pool(graph, protocol="microSPLiT", source=source)
    second = S.ligate_split_pool_round(graph, first, ligase, n=2,
                                       protocol="microSPLiT", source=source)
    third = S.ligate_split_pool_round(graph, second, ligase, n=3,
                                      protocol="microSPLiT", source=source)
    captured = S.capture_split_pool_cdna(graph, third, protocol="microSPLiT")
    amplified = S.switch_and_amplify_split_pool(graph, captured, rt,
                                                protocol="microSPLiT", source=source)

    fragmentation = graph.input("fragmentation_mix", "WGS enzymatic fragmentation mix",
                                MaterialKind.CATALYST)
    fragmented = graph.product("fragmented_cdna", "enzymatically fragmented barcoded cDNA",
                               MaterialKind.MOLECULAR)
    graph.react("fragment", "Fragment amplified cDNA", StepKind.FRAGMENTATION,
                ((amplified, InputRole.TEMPLATE),
                 (fragmentation, InputRole.CATALYST)), (fragmented,))
    prepared = graph.product("prepared_fragments", "end-repaired, dA-tailed cDNA fragments",
                             MaterialKind.MOLECULAR, inferred=True,
                             source=source + " (dA tail inferred from the disclosed 3′-T adapter)")
    graph.react("end_prepare", "Prepare fragment ends for the 3-prime-T adapter",
                StepKind.END_REPAIR,
                ((fragmented, InputRole.TEMPLATE),
                 (fragmentation, InputRole.CATALYST)), (prepared,),
                source=source + " (end repair/dA tail inferred from adapter compatibility)")

    adapter = graph.input("library_adapter", "pre-annealed BC_0243/BC_0244 adapter",
                          MaterialKind.ADAPTER)
    library_ligase = graph.input("library_ligase", "WGS DNA ligase", MaterialKind.CATALYST)
    template = graph.product("adapter_ligated", "adapter-ligated barcoded cDNA fragments",
                             MaterialKind.MOLECULAR)
    graph.react("adapter_ligation", "Ligate the 3-prime-T library adapter", StepKind.LIGATION,
                ((prepared, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (library_ligase, InputRole.CATALYST)), (template,),
                geometry=reported_adapter_end(
                    "BC_0243/BC_0244 adapter", kind="cohesive", polarity="3-prime",
                    terminal_sequence="T",
                    statement="The microSPLiT oligo table and fragmentation method disclose a duplex adapter with a terminal 3′ T."))

    forward = graph.input("pcr_forward", "library PCR forward primer", MaterialKind.OLIGO,
                          state=Construct([Segment("BC_0027", S.MICRO_LIBRARY_FWD)]),
                          source=source + " Supplementary Table S3")
    reverse = graph.input("pcr_reverse", "indexed library PCR reverse primer",
                          MaterialKind.OLIGO,
                          state=Construct([Segment("BC_0076 / TSBC07", S.MICRO_LIBRARY_REV)]),
                          source=source + " Supplementary Table S3")
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase", MaterialKind.CATALYST)
    final = S.final_boundary("microSPLiT")
    library = graph.product("library", final.name, MaterialKind.LIBRARY, state=final,
                            source=source + " Supplementary Table S3")
    graph.react("library_pcr", "Build the indexed sequencing library", StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.exact(
                    final, PCRPrimer("BC_0027", S.MICRO_LIBRARY_FWD),
                    PCRPrimer("BC_0076 / TSBC07", S.MICRO_LIBRARY_REV),
                    statement="Supplementary Table S3 prints the complete terminal PCR pair."))
    graph.sequence(library, SequencingPlan(
        library=final, primers=S.MICRO_SEQ_PRIMERS,
        acquired_roles=("Read 1", "Read 2", "Index 1 (i7)"),
        read_lengths={"Read 1": 74, "Read 2": 86, "Index 1 (i7)": 6},
        read_uses={"Read 1": "include", "Read 2": "extract", "Index 1 (i7)": "demultiplex"}))
    graph.declare_reviewed_scope((
        "fixed_cells", "round1_primer", "reverse_transcriptase", "reverse_transcribe",
        "dna_ligase", "round2_adapter", "round2_splint", "ligate_round2",
        "round3_adapter", "round3_splint", "ligate_round3", "streptavidin", "capture",
        "tso", "template_switch", "preamp_forward", "preamp_reverse",
        "preamp_polymerase", "cdna_pcr", "fragmentation_mix", "fragment", "end_prepare",
        "library_adapter", "library_ligase", "adapter_ligation", "pcr_forward",
        "pcr_reverse", "pcr_polymerase", "library_pcr", "library"), source)
    return graph.require_publishable()
