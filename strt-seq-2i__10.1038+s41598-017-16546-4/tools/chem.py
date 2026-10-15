"""STRT-seq-2i source-owned reaction graph."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "strt-seq__10.1101+gr.110882.110" / "tools"), str(HERE.parents[1] / "lib")]
import strt as S
import illumina as il
from chemdraw import Construct
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind)
from page_parts import _front, _pool, _tagment, seg


def protocol_graph() -> ProtocolGraph:
    source = "https://doi.org/10.1038/s41598-017-16546-4"
    graph = ProtocolGraph("STRT-seq-2i", authority=source)
    switched, reviewed = _front(graph, S.TWO_I)
    indexed = graph.input("indexed_p1b_primer", "DI-P1A-index-P1B primer",
                          MaterialKind.OLIGO,
                          state=Construct([seg("DI-P1A-index-P1B",
                                               S.TWO_I_WELL_PRIMER,
                                               placeholder=True)]))
    polymerase = graph.input("indexed_pcr_polymerase", "KAPA HiFi polymerase",
                             MaterialKind.CATALYST)
    indexed_state = S.two_i_preamp_product()
    grafted = graph.product("well_indexed_cdna", "well-index-grafted cDNA",
                            MaterialKind.MOLECULAR, state=indexed_state)
    graph.react("graft_well_index", "Graft P5 and the well index at the TSO end",
                StepKind.EXTENSION,
                ((switched, InputRole.TEMPLATE),
                 (indexed, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (grafted,))
    p1a_primer = PCRPrimer("DI-PCR-P1A", il.P5[:20])
    p1a = graph.input("p1a_primer", "DI-PCR-P1A", MaterialKind.OLIGO,
                      state=Construct([seg("DI-PCR-P1A", il.P5[:20])]))
    amplified = graph.product("amplified_cdna", "well-indexed amplified cDNA",
                              MaterialKind.MOLECULAR, state=indexed_state)
    graph.react("indexed_pcr", "Amplify the well-indexed cDNA with P1A",
                StepKind.PCR,
                ((grafted, InputRole.TEMPLATE), (p1a, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.symmetric(
                    indexed_state, p1a_primer,
                    statement="After index grafting, P1A has inward sites at both ends."))
    pooled, pool_reviewed = _pool(
        graph, amplified, indexed_state, identifier="subarray_pool",
        label="other well-indexed cDNA from the subarray")
    reviewed.extend(pool_reviewed)
    tagged = _tagment(graph, pooled, S.two_i_library(), reviewed)
    capture = graph.input("streptavidin_beads", "streptavidin beads",
                          MaterialKind.SELECTION)
    captured = graph.product("captured_cdna", "captured biotinylated tagmented cDNA",
                             MaterialKind.MOLECULAR, state=S.two_i_library())
    capture_waste = graph.product("capture_waste", "unbound tagmented cDNA",
                                  MaterialKind.WASTE)
    graph.react("capture", "Capture biotinylated tagmented cDNA",
                StepKind.AFFINITY_CAPTURE,
                ((tagged, InputRole.TEMPLATE),
                 (capture, InputRole.SELECTION)), (captured, capture_waste))
    graph.discard(capture_waste, "not retained by streptavidin capture")
    exosap = graph.input("exosap", "ExoSAP-IT", MaterialKind.CATALYST)
    cleaned = graph.product("adapter_cleaned_cdna", "adapter-cleaned captured cDNA",
                            MaterialKind.MOLECULAR, state=S.two_i_library())
    cleanup_waste = graph.product("adapter_cleanup_waste",
                                  "digested residual adapters", MaterialKind.WASTE)
    graph.react("remove_adapters", "Remove leftover adapters with ExoSAP-IT",
                StepKind.CLEANUP,
                ((captured, InputRole.TEMPLATE),
                 (exosap, InputRole.CATALYST)), (cleaned, cleanup_waste))
    graph.discard(cleanup_waste, "residual adapters digested before elution")
    heat = graph.input("heat_elution", "70 °C strand-elution condition",
                       MaterialKind.CONDITION)
    ss_template = graph.product("single_strand_template",
                                "eluted single-stranded library template",
                                MaterialKind.MOLECULAR, state=S.two_i_library(),
                                inferred=True,
                                source=source + " (P2/barcode junction unresolved)")
    strand_waste = graph.product("captured_strand_waste", "bead-retained strand",
                                 MaterialKind.WASTE)
    graph.react("strand_elution", "Elute the single-stranded library template",
                StepKind.STRAND_SEPARATION,
                ((cleaned, InputRole.TEMPLATE),
                 (heat, InputRole.CONDITION)), (ss_template, strand_waste))
    graph.discard(strand_waste, "remains bead-associated after heat elution")
    final = S.two_i_library()
    p1_primer = PCRPrimer("P1 second-PCR", il.P5[:24])
    p2_primer = PCRPrimer("P2 second-PCR", il.P7)
    p1 = graph.input("p1_second_pcr", p1_primer.name, MaterialKind.OLIGO,
                     state=Construct([seg(p1_primer.name, p1_primer.sequence)]))
    p2 = graph.input("p2_second_pcr", p2_primer.name, MaterialKind.OLIGO,
                     state=Construct([seg(p2_primer.name, p2_primer.sequence)]))
    second_polymerase = graph.input("second_pcr_polymerase", "KAPA HiFi polymerase",
                                   MaterialKind.CATALYST)
    unselected_library = graph.product(
        "unselected_library", "amplified STRT-seq-2i library",
        MaterialKind.MOLECULAR, state=final, inferred=True,
        source=source + " (P2/barcode junction unresolved)")
    graph.react("second_pcr", "Amplify the selected strand with P1/P2",
                StepKind.PCR,
                ((ss_template, InputRole.TEMPLATE), (p1, InputRole.STRUCTURAL),
                 (p2, InputRole.STRUCTURAL),
                 (second_polymerase, InputRole.CATALYST)), (unselected_library,),
                geometry=PCRGeometry.exact(
                    final, p1_primer, p2_primer,
                    statement="Table 1 prints both second-PCR primers."))
    beads = graph.input("size_selection_beads", "AMPure size-selection beads",
                        MaterialKind.SELECTION)
    library = graph.product("library", "size-selected STRT-seq-2i library",
                            MaterialKind.LIBRARY, state=final, inferred=True,
                            source=source + " (P2/barcode junction unresolved)")
    size_waste = graph.product("size_waste", "rejected library sizes",
                               MaterialKind.WASTE)
    graph.react("size_select", "Double-size-select the final library",
                StepKind.SELECTION,
                ((unselected_library, InputRole.TEMPLATE),
                 (beads, InputRole.SELECTION)), (library, size_waste))
    graph.discard(size_waste, "outside the retained AMPure size interval")
    graph.sequence(library, SequencingPlan(
        library=final, primers=S.TWO_I_INDEX_PRIMERS,
        acquired_roles=tuple(p.role for p in S.TWO_I_INDEX_PRIMERS)))
    reviewed.extend([
        "indexed_p1b_primer", "indexed_pcr_polymerase", "graft_well_index",
        "p1a_primer", "indexed_pcr", "streptavidin_beads", "capture", "exosap",
        "remove_adapters", "heat_elution", "strand_elution", "p1_second_pcr",
        "p2_second_pcr", "second_pcr_polymerase", "second_pcr",
        "size_selection_beads", "size_select", "library",
    ])


    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
