"""Source-backed page specifications for targeted 3C and Repli-seq methods."""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import MolecularState, Scene, Segment, Workflow
from chromosome_capture import region

def role(name: str, *, tag=None, molecular_feature=None) -> Segment:
    return Segment(name, f"[{name}]", tag, placeholder=True,
                   feature=molecular_feature)

def molecule(label: str, names: tuple[str, ...], *,
             junctions: tuple[tuple[str, str, str], ...] = ()) -> Scene:
    parts = [region(name) for name in names]
    sc = Scene(); sc.strand("molecule", parts, label=label)
    for left, right, text in junctions:
        sc.junction("molecule", left, right, text)
    sc.labels("molecule")
    return sc

def chain(initial: Scene, steps: tuple[tuple[str, Scene, str], ...]) -> Workflow:
    w = Workflow(MolecularState("Starting material", tuple(initial.rows())))
    for action, state, note in steps:
        w.react(action, state.rows(), name=action, note=note)
    return w


@dataclass(frozen=True)
class Definition:
    title: str
    doi: str
    citation: str
    summary: str
    final_parts: tuple[Segment, ...]
    sequencing: str
    standard_truseq: bool
    caveat: str = ""
