"""Small shared constructors for role-level protocol schematics.

These are used when a source establishes molecular roles and transformations but does
not publish one immutable base-level adapter geometry.  They keep placeholder handling,
junction marking and workflow continuity in one place.
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
             junctions: tuple[tuple[str, str, str], ...] = ()) -> Scene:
    sc = Scene()
    sc.strand("molecule", list(parts), label=label)
    for left, right, kind in junctions:
        sc.junction("molecule", left, right, kind)
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


def role_library(name: str, payload: tuple[Segment, ...], *,
                 index: Segment | None = None) -> Construct:
    parts = [role("Illumina-compatible left end")]
    if index is not None:
        parts.append(index)
    parts.extend(payload)
    parts.append(role("Illumina-compatible right end"))
    return Construct(parts, name=name)
