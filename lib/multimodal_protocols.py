"""Metadata compatibility for locally owned multimodal protocol declarations."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Definition:
    title: str; doi: str; citation: str; summary: str; endpoint: str; notes: str
