"""Reusable TruSeq/Nextera molecular states and bounded workflow operations."""
from __future__ import annotations

from dataclasses import dataclass

import nextera as nx
import illumina as il
from batch_ngs import (nextera_library, seg, truseq_full_length_adapter, truseq_library,
                       truseq_pcr_free_library, truseq_pcr_free_rows)
from chemdraw import (Construct, MolecularState, Scene, Segment, Workflow,
                      complement_segments)
from endprep import (ThreePrimeOverhangAdapter, cohesive_ligation_geometry,
                     dA_tailed_scene, ligation_product_scene, repair_and_dA_tail)
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, SequencingPlan, StepKind

TRUSEQ_Y_END = ThreePrimeOverhangAdapter(
    "TruSeq Y adapter", il.STEM + "T", il.STEM_COMPLEMENT)

def _duplex(name: str, *, bases: str = "N" * 36) -> Scene:
    return Scene.duplex([Segment(name, bases, placeholder=True)], label=name)

def _single(name: str, parts: list[Segment], *, label: str) -> Scene:
    sc = Scene()
    sc.strand(name, parts, label=label)
    sc.labels(name)
    return sc

def _workflow(initial_name: str, initial: Scene,
              steps: tuple[tuple[str, str, Scene, str], ...]) -> Workflow:
    w = Workflow(MolecularState(initial_name, tuple(initial.rows())))
    for action, name, state, note in steps:
        w.react(action, state.rows(), name=name, note=note)
    return w

def _truseq_final(insert_name: str, title: str, *, dual_index: bool = True):
    return truseq_library([seg(insert_name, "N" * 42, placeholder=True)], title,
                          dual_index=dual_index)

def _truseq_pcr_free_final(insert_name: str, title: str, *, adapter):
    return truseq_pcr_free_library([seg(insert_name, "N" * 42, placeholder=True)], title,
                                   adapter=adapter)

def _nextera_final(insert_name: str, title: str):
    return nextera_library([seg(insert_name, "N" * 42, placeholder=True)], title)

def truseq_final_rows(library: Construct) -> list:
    insert = next(s.name for s in library if "insert" in s.name or "cDNA" in s.name)
    return ligation_product_scene(
        library, insert, label=library.name, adapter_group="TruSeq adapter",
        highlight_components=False).rows()

def nextera_final_rows(library: Construct) -> list:
    sc = Scene.duplex(list(library), label=library.name)
    insert = next(s.name for s in library if "insert" in s.name or "DNA" in s.name)
    sc.junction("top", "left mosaic end", insert, "transposition boundary")
    sc.junction("top", insert, "right mosaic end reverse complement",
                "transposition boundary")
    return sc.rows()

def _tagged_intermediate(label: str, *, gaps: bool = False) -> Scene:
    """Heterologous S5/S7 Tn5 product, drawn as the caller's source supports.

    ``gaps=True`` draws the transposition product before gap filling: each transferred
    adapter (entry point + mosaic end) is joined to a genomic 5' end and Tn5's 9-bp
    stagger leaves the nine genomic bases beside each join unpaired.  The
    non-transferred mosaic-end strand is not drawn, because vendor documents do not
    print it.  The legacy default (``gaps=False``) is a closed-duplex simplification
    that omits the gaps; a page using it must say so.
    """
    if not gaps:
        parts = [Segment("S5", nx.S5, "s5"), Segment("left mosaic end", nx.ME, "me"),
                 Segment(label, "N" * 36, placeholder=True),
                 Segment("right mosaic end reverse complement", nx.ME_RC, "me"),
                 Segment("S7 reverse complement", nx.S7_RC, "s7")]
        sc = Scene.duplex(parts, label="tagmented DNA")
        sc.junction("top", "left mosaic end", label, "transposition boundary")
        sc.junction("top", label, "right mosaic end reverse complement",
                    "transposition boundary")
        return sc
    gap = nx.TAGMENTATION_GAP
    core = Segment(label, "N" * 27, placeholder=True)
    left_gap = Segment("left 9-nt end", "N" * gap, placeholder=True)
    right_gap = Segment("right 9-nt end", "N" * gap, placeholder=True)
    top = [Segment("S5", nx.S5, "s5"), Segment("left mosaic end", nx.ME, "me"),
           left_gap, core]
    bottom = [Segment("S7", nx.S7, "s7"), Segment("right mosaic end", nx.ME, "me"),
              right_gap, *complement_segments([core])]
    sc = Scene()
    sc.strand("top", top, label="tagmented DNA")
    sc.anneal("bottom", bottom, to="top", pair=(label + "'", label), label="tagmented DNA")
    sc.junction("top", "left mosaic end", left_gap.name, "transposition boundary")
    sc.junction("bottom", "right mosaic end", right_gap.name, "transposition boundary")
    sc.mark("top", left_gap.name, "unpaired: 9-nt gap opposite")
    sc.mark("bottom", right_gap.name, "unpaired: 9-nt gap opposite")
    sc.labels("top")
    return sc

@dataclass(frozen=True)
class Definition:
    title: str
    notes: str
    source: str
    summary: str
    caveat: str
    final_caption: str
    sequencing_intro: str
    final_kind: str
    insert_name: str
    workflow_builder: object

def _sequence(graph, material, library, primers, endpoint):
    graph.sequence(material, SequencingPlan(
        library=library, primers=tuple(primers), limitation=""))

def _finish_pcr(graph, prefix, template, library, primers, endpoint, reviewed, *, label,
                step_label=None, forward_label=None, reverse_label=None,
                polymerase_label=None, forward_state=None, reverse_state=None,
                primer_source=None, geometry=None, unavailable_reason=None,
                cleanup=None):
    """Legacy vendor index PCR; callers supply printed primers, geometry and cleanup.

    ``geometry`` takes an exact :class:`PCRGeometry` for a kit revision whose primers
    are printed.  ``cleanup`` is an optional ``(step_label, reagent_label,
    discarded_label)`` tuple for a post-PCR bead cleanup before sequencing.
    """
    forward = graph.input(prefix + "_forward", forward_label or label + " index PCR forward primer",
                          MaterialKind.OLIGO, state=forward_state, source=primer_source)
    reverse = graph.input(prefix + "_reverse", reverse_label or label + " index PCR reverse primer",
                          MaterialKind.OLIGO, state=reverse_state, source=primer_source)
    polymerase = graph.input(prefix + "_polymerase", polymerase_label or label + " PCR polymerase",
                             MaterialKind.CATALYST)
    amplified_id = prefix + ("_amplified" if cleanup else "_library")
    product = graph.product(amplified_id, label + (" amplified library" if cleanup
                                                    else " sequencing library"),
                            MaterialKind.MOLECULAR if cleanup else MaterialKind.LIBRARY,
                            state=library)
    graph.react(prefix + "_pcr", step_label or "Complete and index the " + label + " library",
                StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (product,), geometry=geometry or PCRGeometry.unavailable(
                    library,
                    reason=unavailable_reason or "The vendor workflow supports several index-component revisions; the canonical final library and sequencing sites are exact, while one universal PCR-primer pair is not."))
    reviewed.extend((prefix + "_forward", prefix + "_reverse", prefix + "_polymerase",
                     prefix + "_pcr", amplified_id))
    if cleanup:
        cleanup_label, reagent_label, discarded_label = cleanup
        beads = graph.input(prefix + "_cleanup_beads", reagent_label, MaterialKind.SELECTION)
        cleaned = graph.product(prefix + "_library", label + " sequencing library",
                                MaterialKind.LIBRARY, state=library)
        discarded = graph.product(prefix + "_cleanup_waste", discarded_label,
                                  MaterialKind.WASTE)
        graph.react(prefix + "_cleanup", cleanup_label, StepKind.CLEANUP,
                    ((product, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                    (cleaned, discarded), changes_molecule=False)
        graph.discard(discarded, "left in the bead supernatant")
        reviewed.extend((prefix + "_cleanup_beads", prefix + "_cleanup", prefix + "_library"))
        product = cleaned
    _sequence(graph, product, library, primers, endpoint)
    return product
