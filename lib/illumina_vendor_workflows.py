"""Reusable TruSeq/Nextera molecular states and bounded workflow operations."""
from __future__ import annotations

from dataclasses import dataclass

import nextera as nx
import illumina as il
from batch_ngs import (nextera_library, seg, truseq_full_length_adapter, truseq_library,
                       truseq_pcr_free_library, truseq_pcr_free_rows)
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow
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

def _tagged_intermediate(label: str) -> Scene:
    parts = [Segment("S5", nx.S5, "s5"), Segment("left mosaic end", nx.ME, "me"),
             Segment(label, "N" * 36, placeholder=True),
             Segment("right mosaic end reverse complement", nx.ME_RC, "me"),
             Segment("S7 reverse complement", nx.S7_RC, "s7")]
    sc = Scene.duplex(parts, label="tagmented DNA")
    sc.junction("top", "left mosaic end", label, "transposition boundary")
    sc.junction("top", label, "right mosaic end reverse complement",
                "transposition boundary")
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

def _finish_pcr(graph, prefix, template, library, primers, endpoint, reviewed, *, label):
    forward = graph.input(prefix + "_forward", label + " index PCR forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input(prefix + "_reverse", label + " index PCR reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input(prefix + "_polymerase", label + " PCR polymerase",
                             MaterialKind.CATALYST)
    product = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=library)
    graph.react(prefix + "_pcr", "Complete and index the " + label + " library",
                StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (product,), geometry=PCRGeometry.unavailable(
                    library,
                    reason="The vendor workflow supports several index-component revisions; the canonical final library and sequencing sites are exact, while one universal PCR-primer pair is not."))
    _sequence(graph, product, library, primers, endpoint)
    reviewed.extend((prefix + "_forward", prefix + "_reverse", prefix + "_polymerase",
                     prefix + "_pcr", prefix + "_library"))
    return product
