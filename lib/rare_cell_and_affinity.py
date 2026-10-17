"""Bounded scene and workflow helpers for selection assays."""
from __future__ import annotations

from chemdraw import MolecularState, Scene, Segment, Workflow


def _scene(label: str, *parts: Segment) -> Scene:
    sc = Scene()
    sc.strand("molecule", list(parts), label=label)
    sc.labels("molecule")
    return sc

def _workflow(initial: Scene, steps: tuple[tuple[str, Scene, str], ...]) -> Workflow:
    w = Workflow(MolecularState("Starting material", tuple(initial.rows())))
    for action, state, note in steps:
        w.react(action, state.rows(), name=action, note=note)
    return w
