"""Canonical paper display metadata for catalogue identifiers needing repair.

Most entries predate strict DOI resolution and encode a malformed or noncanonical
identifier. A current preprint may also temporarily be absent from Crossref and NCBI.
Keeping the correction here prevents the public index from losing an otherwise verified
paper citation.
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
    "chip-seq__10.1038+nmeth1068": PaperReference(
        "Robertson", "2007", "Nature Methods", "https://doi.org/10.1038/nmeth1068"),
    "proximity-network-assay__10.1101+2025.06.19.660329": PaperReference(
        "Karlsson", "2025", "bioRxiv",
        "https://doi.org/10.1101/2025.06.19.660329"),
}


# Catalogue fields for an ours-only row (a protocol upstream does not list) that no
# resolver supplies. fetch_scg_lib_structs.py fills these into scg_lib_structs.tsv after
# NCBI's ID converter has run; a value here never replaces a resolved one, and a
# resolved value that disagrees with it stops the build. Keyed by protocol directory.
CATALOGUE_FIELDS = {
    "chip-seq__10.1038+nmeth1068": {"citation_text": "ChIP-seq"},
}
