"""Reusable base-state transformations used by conversion and damage protocols.

These functions model the chemically informative base, rather than pretending that an
entire library sequence is known.  Unsupported inputs fail at construction time.
"""
from __future__ import annotations

from dataclasses import dataclass


CYTOSINES = {"C", "5mC", "5hmC", "5caC", "DHU", "U", "T", "s4U", "IAA-s4U"}


@dataclass(frozen=True)
class BasePath:
    states: tuple[str, ...]

    def __post_init__(self) -> None:
        if len(self.states) < 2 or any(x not in CYTOSINES for x in self.states):
            raise ValueError("a base-conversion path needs at least two supported states")

    def text(self) -> str:
        return " → ".join(self.states)


def taps_path(base: str) -> BasePath:
    """TAPS chemistry and PCR readout for one cytosine state."""
    if base in {"5mC", "5hmC"}:
        return BasePath((base, "5caC", "DHU", "T"))
    if base == "C":
        return BasePath(("C", "C"))
    raise ValueError("TAPS input must be C, 5mC or 5hmC")


def bisulfite_path(*, protected: bool) -> BasePath:
    return BasePath(("5mC", "C")) if protected else BasePath(("C", "U", "T"))


def slam_path(*, labelled: bool) -> BasePath:
    """SLAM-seq base path through IAA derivatization and RT readout."""
    return BasePath(("s4U", "IAA-s4U", "C")) if labelled else BasePath(("U", "T"))


@dataclass(frozen=True)
class PartialUDGProduct:
    """Surviving strand after USER treatment, before end repair.

    Terminal uracils are retained by protocol design; internal uracils are cleavage
    sites and cannot remain in a constructed product.
    """
    left_terminal: str
    interior: str
    right_terminal: str

    def __post_init__(self) -> None:
        if self.left_terminal not in {"U", "C"} or self.right_terminal not in {"U", "C"}:
            raise ValueError("terminal damage state must be U or C")
        if "U" in self.interior:
            raise ValueError("partial-UDG product cannot retain an internal uracil")

    def text(self) -> str:
        return self.left_terminal + self.interior + self.right_terminal
