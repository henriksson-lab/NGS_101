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
from protocol_graph import (Disposition, InputRole, MaterialKind, ProtocolGraph,
                            SequencingPlan, SourceClaim, SourceDiscrepancy, StepKind,
                            Terminal)
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
    capture = graph.input("streptavidin_beads", "Dynabeads MyOne Streptavidin C1 in BB buffer",
                          MaterialKind.SELECTION)
    captured = graph.product("captured_cdna", "captured biotinylated tagmented cDNA",
                             MaterialKind.MOLECULAR, state=S.two_i_library())
    capture_waste = graph.product("capture_waste",
                                  "unbound tagmented fragments without a biotinylated end",
                                  MaterialKind.WASTE)
    # Table 1 prints 5'Bio on STRT-P1-T31 and DI-PCR-P1A as well as on the
    # TSO-end primers, so oligo-dT-end terminal fragments are also candidates
    # for streptavidin capture.  No sequence state is claimed for them.
    dt_end = graph.product("captured_dt_end_fragments",
                           "candidate co-captured biotinylated oligo-dT-end fragments",
                           MaterialKind.MOLECULAR, inferred=True,
                           source=source + " Table 1 (5′Bio on STRT-P1-T31 and DI-PCR-P1A)")
    graph.react("capture", "Capture biotinylated tagmented cDNA (RT 15 min)",
                StepKind.AFFINITY_CAPTURE,
                ((tagged, InputRole.TEMPLATE),
                 (capture, InputRole.SELECTION)), (captured, dt_end, capture_waste))
    graph.discard(capture_waste, "fragments without a 5′-biotinylated end are not retained "
                                 "by streptavidin capture and are removed in the TNT washes")
    graph.mark_terminal(Terminal(dt_end.id, Disposition.UNRESOLVED,
        "Candidate only: both cDNA ends carry 5′ biotin in Table 1, and the inspected "
        "Methods describe no separation of oligo-dT-end fragments at capture. P1_2nd_PCR "
        "differs from the oligo-dT-end AATGATACGGCGACCACCGATCG sequence in its last four "
        "3′ bases, and analysis discards reads ending in poly(A) with <25 remaining bases; "
        "neither establishes physical removal or carryover."))
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
    bind_beads = graph.input("ampure_binding_beads", "1.5X AMPure beads",
                             MaterialKind.SELECTION)
    bound_template = graph.product("bead_bound_template",
                                   "AMPure-bound single-stranded library template",
                                   MaterialKind.MOLECULAR, state=S.two_i_library(),
                                   inferred=True,
                                   source=source + " (P2/barcode junction unresolved)")
    bind_waste = graph.product("ampure_binding_supernatant",
                               "unbound supernatant", MaterialKind.WASTE)
    graph.react("ampure_bind", "Bind the eluate to 1.5X AMPure beads; PCR on beads",
                StepKind.CLEANUP,
                ((ss_template, InputRole.TEMPLATE),
                 (bind_beads, InputRole.SELECTION)), (bound_template, bind_waste))
    graph.discard(bind_waste, "supernatant removed before resuspending beads in second-PCR mix")
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
    graph.react("second_pcr", "Amplify the selected strand with P1/P2 (8 cycles)",
                StepKind.PCR,
                ((bound_template, InputRole.TEMPLATE), (p1, InputRole.STRUCTURAL),
                 (p2, InputRole.STRUCTURAL),
                 (second_polymerase, InputRole.CATALYST)), (unselected_library,),
                geometry=PCRGeometry.exact(
                    final, p1_primer, p2_primer,
                    statement="Table 1 prints both second-PCR primers."))
    cleanup_beads = graph.input("pcr_cleanup_beads", "0.7X AMPure beads",
                                MaterialKind.SELECTION)
    cleaned_library = graph.product(
        "pcr_cleaned_library", "0.7X-cleaned STRT-seq-2i library",
        MaterialKind.MOLECULAR, state=final, inferred=True,
        source=source + " (P2/barcode junction unresolved)")
    pcr_cleanup_waste = graph.product("pcr_cleanup_waste",
                                      "unbound primers and short products",
                                      MaterialKind.WASTE)
    graph.react("pcr_cleanup", "Clean the transferred PCR supernatant with 0.7X AMPure",
                StepKind.CLEANUP,
                ((unselected_library, InputRole.TEMPLATE),
                 (cleanup_beads, InputRole.SELECTION)), (cleaned_library, pcr_cleanup_waste))
    graph.discard(pcr_cleanup_waste, "not retained by 0.7X AMPure binding")
    beads = graph.input("size_selection_beads", "0.5X then 1X AMPure size-selection beads",
                        MaterialKind.SELECTION)
    library = graph.product("library", "size-selected STRT-seq-2i library",
                            MaterialKind.LIBRARY, state=final, inferred=True,
                            source=source + " (P2/barcode junction unresolved)")
    size_waste = graph.product("size_waste", "rejected library sizes",
                               MaterialKind.WASTE)
    graph.react("size_select", "Double-size-select the final library (0.5X, then 1X)",
                StepKind.SELECTION,
                ((cleaned_library, InputRole.TEMPLATE),
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
        "captured_dt_end_fragments", "ampure_binding_beads", "ampure_bind",
        "pcr_cleanup_beads", "pcr_cleanup",
    ])
    graph.annotate(SourceDiscrepancy(
        "p1_second_pcr", "Second-PCR primer names and concentration",
        (SourceClaim("Sci Rep 2017 Table 1, PCR 2",
                     "P1_2nd_PCR AATGATACGGCGACCACCGAGATC and P2_2nd_PCR CAAGCAGAAGACGGCATACGAGAT."),
         SourceClaim("Sci Rep 2017 Methods, 'Tagmentation and isolation of 5′ fragments'",
                     "Second PCR mix with '200 mM 4K-P1_2ndPCR, 200 nM P2_4K_2ndPCR'.")),
        "The model uses the Table 1 sequences for both second-PCR primers and does not draw a concentration.",
        "Table 1 is the paper's only printed sequence for each primer, and the Methods state that primer sequences for these steps are given in Table 1.",
        "The paper does not define the '4K' names; that they denote the Table 1 oligos and that '200 mM' means 200 nM are readings, not printed facts."))


    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
