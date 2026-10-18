"""Constructors for declared molecular roles and legacy schematic shortcuts.

An unassigned role records the model's scope; it does not establish source absence.
Supplied parts can share placeholder handling, drawing and workflow continuity.
Legacy ``role_library`` and ``da_tailed`` choose platform or end representations;
their caller reviews are tracked in catalogue/chemistry_decoupling.tsv. Do not add
callers to those shortcuts while their source choices are being migrated.
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


def da_tailed(label: str, insert: str = "end-repaired DNA insert") -> Scene:
    """Role-level product after dA-tailing and before adapter ligation."""
    return molecule(label, (role("3-prime dA"), role(insert), role("3-prime dA")))


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
