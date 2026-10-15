"""Reusable state, scene, and library primitives for RNA and cell assays."""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry, reported_adapter_end
from batch_ngs import nextera_library, seg
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature, revcomp
from circular import CircularizationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            SourceCoverage, SourceItem, StepKind)

def token(name: str, *, tag: str | None = None, length: int = 18,
          molecular_feature=None, inferred: bool = False) -> Segment:
    return Segment(name, "X" * length, tag, placeholder=True,
                   feature=molecular_feature, inferred=inferred)

def unknown_region(name: str, *, tag: str | None = None, molecular_feature=None) -> Segment:
    """A source-established role whose exact sequence or length is not established."""
    return Segment(name, f"[{name}]", tag, placeholder=True,
                   feature=molecular_feature)

def single(label: str, parts: tuple[Segment, ...], *, note: str = "",
           junctions: tuple[tuple[str, str, str], ...] = ()) -> Scene:
    """One selectable strand with checked, named covalent junctions."""
    sc = Scene()
    sc.strand("molecule", list(parts), label=label)
    for left, right, caption in junctions:
        sc.junction("molecule", left, right, caption)
    if note:
        sc.note("molecule", note)
    sc.labels("molecule")
    return sc

def branches(label: str, rows: tuple[tuple[str, tuple[Segment, ...]], ...]) -> Scene:
    sc = Scene()
    for row, parts in rows:
        sc.strand(row, list(parts), label=label)
        sc.labels(row)
    return sc

def chain(initial_name: str, initial: Scene,
          transitions: tuple[tuple[str, Scene, str], ...]) -> Workflow:
    w = Workflow(MolecularState(initial_name, tuple(initial.rows())))
    for action, state, note in transitions:
        w.react(action, state.rows(), name=action, note=note)
    return w

def custom_library(name: str, payload: tuple[Segment, ...]) -> Construct:
    return Construct([token("Illumina-compatible left arm"), *payload,
                      token("Illumina-compatible right arm")], name=name)

@dataclass(frozen=True)
class Definition:
    title: str
    doi: str
    citation: str
    notes: str
    summary: str
    endpoint: str
    caveat: str = ""
