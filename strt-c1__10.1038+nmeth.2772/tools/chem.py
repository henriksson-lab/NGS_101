"""STRT-C1 source-owned reaction graph."""
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
    source = "https://doi.org/10.1038/nmeth.2772"
    graph = ProtocolGraph("STRT-C1", authority=source)
    switched, reviewed = _front(graph, S.C1)
    primer = PCRPrimer("C1-P1-PCR-2", S.C1_PCR)
    preamp_oligo = graph.input("preamp_primer", "C1-P1-PCR-2",
                               MaterialKind.OLIGO,
                               state=Construct([seg(primer.name, primer.sequence)]))
    polymerase = graph.input("preamp_polymerase", "cDNA-PCR polymerase",
                             MaterialKind.CATALYST)
    preamp_state = S.c1_preamp_product()
    amplified = graph.product("amplified_cdna", "amplified C1 cDNA",
                              MaterialKind.MOLECULAR, state=preamp_state)
    graph.react("preamp", "Amplify full-length C1 cDNA", StepKind.PCR,
                ((switched, InputRole.TEMPLATE),
                 (preamp_oligo, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (amplified,),
                geometry=PCRGeometry.symmetric(
                    preamp_state, primer,
                    statement="C1-P1-PCR-2 binds the shared C1 handles."))
    tagged = _tagment(graph, amplified, S.c1_library(), reviewed,
                      inferred=True, source=source)
    pooled, pool_reviewed = _pool(
        graph, tagged, S.c1_library(), identifier="pooled_tagged_cdna",
        label="other chamber-indexed tagmented cDNA")
    reviewed.extend(pool_reviewed)
    capture = graph.input("streptavidin_beads", "streptavidin beads",
                          MaterialKind.SELECTION, inferred=True,
                          source=source + " (selection details unavailable in accessible source)")
    bound = graph.product("captured_cdna", "captured biotinylated cDNA",
                          MaterialKind.MOLECULAR, state=S.c1_library(),
                          inferred=True,
                          source=source + " (selection reconstructed from upstream source)")
    capture_waste = graph.product("capture_waste", "unbound cDNA",
                                  MaterialKind.WASTE)
    graph.react("capture", "Capture biotinylated cDNA", StepKind.AFFINITY_CAPTURE,
                ((pooled, InputRole.TEMPLATE),
                 (capture, InputRole.SELECTION)), (bound, capture_waste))
    graph.discard(capture_waste, "not retained by streptavidin capture")
    pvui = graph.input("pvui", "PvuI", MaterialKind.CATALYST, inferred=True,
                       source=source + " (selection details unavailable in accessible source)")
    selected = graph.product("selected_five_prime_fragments",
                             "selected transcript-5′ fragments",
                             MaterialKind.MOLECULAR, state=S.c1_library(),
                             inferred=True,
                             source=source + " (selection reconstructed from upstream source)")
    digest_waste = graph.product("released_three_prime_fragments",
                                 "released transcript-3′ fragments",
                                 MaterialKind.WASTE)
    graph.react("pvui_select", "Release oligo-dT-end fragments with PvuI",
                StepKind.DIGESTION,
                ((bound, InputRole.TEMPLATE), (pvui, InputRole.CATALYST)),
                (selected, digest_waste))
    graph.discard(digest_waste, "released from the captured material")
    alkali = graph.input("alkali", "alkaline strand-elution condition",
                         MaterialKind.CONDITION, inferred=True,
                         source=source + " (elution details unavailable in accessible source)")
    final = S.c1_library()
    library = graph.product("library", "single-stranded STRT-C1 library",
                            MaterialKind.LIBRARY, state=final, inferred=True,
                            source=source + " (strand selection reconstructed)")
    complement_waste = graph.product("complement_waste", "retained complementary strand",
                                     MaterialKind.WASTE)
    graph.react("strand_elution", "Elute the non-biotinylated library strand",
                StepKind.STRAND_SEPARATION,
                ((selected, InputRole.TEMPLATE),
                 (alkali, InputRole.CONDITION)), (library, complement_waste))
    graph.discard(complement_waste, "remains associated with the captured strand")
    graph.sequence(library, SequencingPlan(
        library=final, primers=S.C1_PRIMERS,
        acquired_roles=tuple(p.role for p in S.C1_PRIMERS)))
    reviewed.extend([
        "preamp_primer", "preamp_polymerase", "preamp", "streptavidin_beads",
        "capture", "pvui", "pvui_select", "alkali", "strand_elution", "library",
    ])

    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
