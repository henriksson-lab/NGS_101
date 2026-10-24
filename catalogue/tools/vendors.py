"""Commercial-protocol vendor names used by the public catalogue.

The mapping is deliberately executable rather than hand-copied into index rows: every
commercial protocol without a defining paper must resolve here, so adding a new kit
cannot silently leave the paper/company column blank.
"""
from __future__ import annotations


EXACT = {
    "topo-ta-cloning": "Thermo Fisher Scientific / Invitrogen",
    "q5-site-directed-mutagenesis": "New England Biolabs",
    "quikchange-ii-site-directed-mutagenesis": "Agilent Technologies",
    "surecell-3-wta-for-ddseq": "Illumina / Bio-Rad",
    "dnbseq-library-circularization": "MGI Tech",
    "repli-g-single-cell": "QIAGEN",
    "picoplex-dna-seq": "Takara Bio",
    "sureplex": "Illumina",
}

PREFIXES = (
    ("10x-", "10x Genomics"),
    ("oxford-nanopore-", "Oxford Nanopore Technologies"),
    ("pacbio-", "PacBio"),
    ("nebnext-", "New England Biolabs"),
    ("sureselect-", "Agilent Technologies"),
    ("illumina-", "Illumina"),
    ("nextera-", "Illumina"),
    ("truseq-", "Illumina"),
)


def company_for(directory: str) -> str:
    """Return the declared vendor, failing closed for an unclassified commercial kit."""
    if directory in EXACT:
        return EXACT[directory]
    for prefix, company in PREFIXES:
        if directory.startswith(prefix):
            return company
    raise ValueError(f"commercial protocol {directory!r} has no vendor classification")
