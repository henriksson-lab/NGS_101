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
    graph.react("preamp", "Amplify full-length C1 cDNA: 21 cycles", StepKind.PCR,
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
    final = S.c1_library()
    library = graph.product("library", "STRT-C1 reference library; strand recovery unresolved", MaterialKind.LIBRARY,
                            state=final, inferred=True, source=source+"; exact post-tagmentation recovery remains unverified")
    graph.react("recover_library", "Recover the STRT-C1 library without library PCR", StepKind.MAGIC,
                ((pooled, InputRole.TEMPLATE),), (library,),
                reason="The published Supplementary Note establishes tagmentation and omission of library PCR, but the inspected source does not establish the previously assumed streptavidin/PvuI/alkaline selection sequence. Main Online Methods remain deferred; an upstream drawing is not evidence.")
    graph.sequence(library, SequencingPlan(library=final, primers=S.C1_PRIMERS,
        limitation="Reference primer sites are conditional; run acquisition and exact strand recovery remain unassigned in this bounded pass."))
    reviewed.extend(["preamp_primer", "preamp_polymerase", "preamp", "recover_library", "library"])

    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()
