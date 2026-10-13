"""Shared typed graph for single-cell nuclease-fragment library protocols."""
from __future__ import annotations

from adapter_ligation import reported_adapter_end
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)


def protocol_graph(m, *, name: str, source: str, nuclease: str,
                   inferred_adapters: bool, final_amplification: bool) -> ProtocolGraph:
    graph = ProtocolGraph(name, authority=source)
    cell = graph.input("cell", "one lysed cell", MaterialKind.BIOLOGICAL, primary=True)
    enzyme = graph.input("nuclease", nuclease, MaterialKind.CATALYST)
    digested = graph.product("digested_fragments", f"{nuclease}-released genomic fragments",
                             MaterialKind.MOLECULAR, state=m.digested_fragment_scene())
    graph.react("digest", f"Digest one lysed cell with {nuclease}", StepKind.DIGESTION,
                ((cell, InputRole.TEMPLATE), (enzyme, InputRole.CATALYST)), (digested,))
    inferred_source = source + " (post-digest chemistry inferred from named adapter reagents)"
    repair = graph.input("end_repair_mix", "end-repair and dA-tailing enzymes",
                         MaterialKind.CATALYST, inferred=inferred_adapters,
                         source=inferred_source if inferred_adapters else source)
    tailed = graph.product("a_tailed_fragments", "end-repaired dA-tailed fragments",
                           MaterialKind.MOLECULAR, state=m.a_tailed_fragment_scene(),
                           inferred=inferred_adapters,
                           source=inferred_source if inferred_adapters else source)
    graph.react("end_prepare", "End-repair fragments and add one 3-prime dA",
                StepKind.END_REPAIR,
                ((digested, InputRole.TEMPLATE), (repair, InputRole.CATALYST)),
                (tailed,), source=inferred_source if inferred_adapters else source)
    adapters = graph.input("adapters", "forked single-index Illumina adapter",
                           MaterialKind.ADAPTER, state=m.adapter_scene(),
                           inferred=inferred_adapters,
                           source=inferred_source if inferred_adapters else source)
    ligase = graph.input("ligase", "DNA ligase", MaterialKind.CATALYST,
                         inferred=inferred_adapters,
                         source=inferred_source if inferred_adapters else source)
    ligated = graph.product("ligated_fragments", "adapter-ligated fragments",
                            MaterialKind.MOLECULAR, inferred=inferred_adapters,
                            source=inferred_source if inferred_adapters else source)
    geometry = reported_adapter_end(
        "forked Illumina adapter", kind="cohesive", polarity="3-prime",
        terminal_sequence="T", substrate_scenes=(m.a_tailed_fragment_scene(),),
        statement=("The named universal adapter is reconstructed with a 3-prime dT matched to the dA-tailed fragment."
                   if inferred_adapters else
                   "The cited Illumina Multiplexing adapter presents a 3-prime dT matched to the dA-tailed fragment."))
    graph.react("adapter_ligation", "Ligate the forked adapter", StepKind.LIGATION,
                ((tailed, InputRole.TEMPLATE), (adapters, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (ligated,), geometry=geometry,
                source=inferred_source if inferred_adapters else source)
    p5 = graph.input("p5_primer", "P5-side PCR primer", MaterialKind.OLIGO,
                     state=tuple(m.p5_primer()), inferred=inferred_adapters,
                     source=inferred_source if inferred_adapters else source)
    p7 = graph.input("indexed_p7_primer", "indexed P7-side PCR primer",
                     MaterialKind.OLIGO, state=tuple(m.indexed_p7_primer()),
                     inferred=inferred_adapters,
                     source=inferred_source if inferred_adapters else source)
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    amplified = graph.product("indexed_amplicons", "indexed amplified fragments",
                              MaterialKind.MOLECULAR, state=m.final_library(),
                              inferred=inferred_adapters,
                              source=inferred_source if inferred_adapters else source)
    graph.react("indexed_pcr", "Index-amplify adapter-ligated fragments", StepKind.PCR,
                ((ligated, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                 (p7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (amplified,), geometry=PCRGeometry.unavailable(
                    m.final_library(), reason="The ligated pre-PCR intermediate is not represented base-for-base."))
    selector = graph.input("size_selection", "library size-selection reagent",
                           MaterialKind.SELECTION)
    selected_kind = MaterialKind.MOLECULAR if final_amplification else MaterialKind.LIBRARY
    selected_id = "selected_amplicons" if final_amplification else "library"
    selected = graph.product(selected_id, "size-selected indexed fragments", selected_kind,
                             state=m.final_library(), inferred=inferred_adapters,
                             source=inferred_source if inferred_adapters else source)
    rejected = graph.product("rejected_sizes", "rejected library-size fractions",
                             MaterialKind.WASTE)
    graph.react("size_select", "Select the desired library-size fraction",
                StepKind.SELECTION,
                ((amplified, InputRole.TEMPLATE), (selector, InputRole.SELECTION)),
                (selected, rejected))
    graph.discard(rejected, "off-size molecules are excluded from sequencing")
    if final_amplification:
        library = graph.product("library", m.final_library().name, MaterialKind.LIBRARY,
                                state=m.final_library())
        graph.react("final_pcr", "Amplify the selected library with P5 and P7 primers",
                    StepKind.PCR,
                    ((selected, InputRole.TEMPLATE), (p5, InputRole.STRUCTURAL),
                     (p7, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                    (library,), geometry=PCRGeometry.unavailable(
                        m.final_library(), reason="The selected molecules already carry the verified final primer sites."))
    else:
        library = selected
    graph.sequence(library, SequencingPlan(
        library=m.final_library(), primers=m.SEQ_PRIMERS,
        acquired_roles=tuple(p.role for p in m.SEQ_PRIMERS),
        read_uses={p.role: ("demultiplex" if "Index" in p.role else "include")
                   for p in m.SEQ_PRIMERS},
        limitation=("Adapter bases are a marked TruSeq reconstruction."
                    if inferred_adapters else "")))
    scope = ["cell", "nuclease", "digest", "end_repair_mix", "end_prepare",
             "adapters", "ligase", "adapter_ligation", "p5_primer",
             "indexed_p7_primer", "pcr_polymerase", "indexed_pcr",
             "size_selection", "size_select"]
    if final_amplification:
        scope.append("final_pcr")
    scope.append("library")
    graph.declare_reviewed_scope(tuple(scope), source)
    return graph.require_publishable()
