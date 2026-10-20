"""Page metadata compatibility for the remaining local declaration."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Definition:
    """Legacy page metadata still declared locally by SPLASH."""
    title: str
    doi: str
    source_label: str
    summary: str
    library: str
    payload: tuple[str, ...]
    sequencing: str
    modality: str
    caveat: str = ""
