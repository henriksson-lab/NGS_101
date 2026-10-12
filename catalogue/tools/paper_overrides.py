"""Canonical paper display metadata for legacy catalogue identifiers needing repair.

These three directories predate strict DOI resolution and encode a malformed or
noncanonical identifier.  Keeping the correction here prevents the public index from
linking to a dead DOI while a future URL migration can rename the directories safely.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PaperReference:
    first_author: str
    year: str
    journal: str
    href: str


BY_DIRECTORY = {
    "chip-seq__10.1038+nmeth.1068": PaperReference(
        "Robertson", "2007", "Nature Methods", "https://doi.org/10.1038/nmeth1068"),
    "pro-seq__10.1101+gr.150027.112": PaperReference(
        "Kwak", "2013", "Science", "https://doi.org/10.1126/science.1229386"),
    "cutac__10.1186+s13059-020-02266-4": PaperReference(
        "Henikoff", "2020", "eLife", "https://doi.org/10.7554/eLife.63274"),
}
