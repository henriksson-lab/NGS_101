"""Shared role-level models for targeted immune, DNA–RNA and sensor readouts."""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

from chemdraw import Construct, Segment, feature
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import chain, molecule, role, role_library

@dataclass(frozen=True)
class Definition:
    title: str
    doi: str
    citation: str
    summary: str
    endpoint: str
    notes: str

def _authority(d: Definition) -> str:
    return f"{d.citation}; doi:{d.doi}"

def _unavailable(final: Construct, reason: str) -> PCRGeometry:
    return PCRGeometry.unavailable(final, reason=reason)
