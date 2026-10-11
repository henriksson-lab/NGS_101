"""Read the manually refreshed OpenAlex citation snapshot.

Citation counts are discovery metadata, not protocol evidence. The site build consumes
only the checked-in snapshot and never contacts OpenAlex.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from datetime import date
from pathlib import Path


TSV = Path(__file__).resolve().parents[1] / "citations.tsv"
COLUMNS = ["doi", "citations", "openalex_id", "retrieved"]


@dataclass(frozen=True)
class Citation:
    doi: str
    citations: int | None
    openalex_id: str
    retrieved: str


def load(path: Path = TSV, expected_dois: set[str] | None = None) -> dict[str, Citation]:
    """Return DOI-keyed records, rejecting malformed or incomplete snapshots.

    When the site supplies its defining DOI set, every DOI must have a row, including an
    explicitly unresolved row.  Adding a protocol therefore cannot silently omit its
    citation metadata: the maintainer is directed to refresh the checked-in snapshot.
    """
    if not path.exists():
        return {}
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames != COLUMNS:
            raise ValueError(f"{path}: columns are {reader.fieldnames}, expected {COLUMNS}")
        rows = list(reader)
    out: dict[str, Citation] = {}
    for line, row in enumerate(rows, 2):
        doi = row["doi"].strip().lower()
        if not doi or doi in out:
            raise ValueError(f"{path}:{line}: missing or duplicate DOI {doi!r}")
        raw = row["citations"].strip()
        if raw and not raw.isdigit():
            raise ValueError(f"{path}:{line}: invalid citation count {raw!r}")
        if bool(raw) != bool(row["openalex_id"].strip()):
            raise ValueError(
                f"{path}:{line}: citation count and OpenAlex ID must both be present or absent"
            )
        retrieved = row["retrieved"].strip()
        try:
            date.fromisoformat(retrieved)
        except ValueError as exc:
            raise ValueError(f"{path}:{line}: invalid retrieval date {retrieved!r}") from exc
        out[doi] = Citation(doi, int(raw) if raw else None,
                            row["openalex_id"].strip(), retrieved)
    if expected_dois is not None:
        expected = {d.strip().lower() for d in expected_dois if d.strip()}
        missing = expected - set(out)
        extra = set(out) - expected
        if missing or extra:
            raise ValueError(
                f"{path}: DOI snapshot mismatch; missing={sorted(missing)}, "
                f"extra={sorted(extra)}. Run python3 tools/fetch_citations.py after "
                "adding or removing a protocol with a defining paper."
            )
    return out
