"""Constructors for explicitly supplied molecular roles and workflow views.

An unassigned role records model scope; it does not establish source absence.
"""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import Construct, MolecularState, Row, Scene, Segment, Workflow, circle_rows


@dataclass(frozen=True)
class Drawing:
    """A checked row drawing that can be used as a workflow molecular state."""
    data: tuple[Row, ...]

    def rows(self) -> list[Row]:
        return list(self.data)


def role(name: str, *, tag: str | None = None, length_bp: int | None = None,
         molecular_feature=None, inferred: bool = False) -> Segment:
    token = f"[{length_bp}-nt {name}]" if length_bp else f"[{name}]"
    return Segment(name, token, tag, placeholder=True, length_bp=length_bp,
                   feature=molecular_feature, inferred=inferred)


def molecule(label: str, parts: tuple[Segment, ...], *,
             junctions: tuple[tuple[str, str, str], ...] = (),
             components: tuple[tuple[str, str, str], ...] = ()) -> Scene:
    """Join ordered parts on one molecular strand.

    A cell, reaction mixture or independent molecular population needs separate
    ``Scene.strand`` rows; putting its components here would join them covalently.
    """
    sc = Scene()
    sc.strand("molecule", list(parts), label=label)
    for left, right, kind in junctions:
        sc.junction("molecule", left, right, kind)
    for first, last, identity in components:
        sc.component("molecule", first, last, identity)
    sc.labels("molecule")
    return sc


def closed_circle(label: str, parts: tuple[Segment, ...], *,
                  closure_label: str = "covalently closed") -> Drawing:
    """A readable circular molecule whose last-to-first closure is model-derived."""
    con = Construct(parts, name=label)
    return Drawing(tuple(circle_rows(con, closure_label)))


def chain(initial: Scene, steps: tuple[tuple[str, Scene, str], ...]) -> Workflow:
    workflow = Workflow(MolecularState("Starting material", tuple(initial.rows())))
    for action, state, note in steps:
        workflow.react(action, state.rows(), name=action, note=note)
    return workflow
