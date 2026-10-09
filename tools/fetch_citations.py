#!/usr/bin/env python3
"""Refresh the static OpenAlex citation snapshot for defining papers.

This command is deliberately manual; normal documentation and site builds never access
the network. Casual keyless use is supported by OpenAlex. Set OPENALEX_API_KEY for a
larger request budget.

    python3 tools/fetch_citations.py
    OPENALEX_API_KEY=... python3 tools/fetch_citations.py
    python3 tools/fetch_citations.py --out /tmp/citations.tsv
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import tempfile
import urllib.error
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "catalogue" / "tools"))

import catalogue  # noqa: E402
from citations import COLUMNS, TSV  # noqa: E402


API = "https://api.openalex.org/works"
BATCH = 50                 # below OpenAlex's documented 100-value OR limit
USER_AGENT = "NGS-101-citation-snapshot/1.0 (https://github.com/henriksson-lab/NGS_101)"


def batches(values: list[str], size: int = BATCH):
    for start in range(0, len(values), size):
        yield values[start:start + size]


def request_batch(dois: list[str], api_key: str = "") -> list[dict]:
    query = {
        "filter": "doi:" + "|".join("https://doi.org/" + d for d in dois),
        "per_page": "100",
        "select": "id,doi,cited_by_count",
    }
    url = API + "?" + urllib.parse.urlencode(query)
    headers = {"User-Agent": USER_AGENT}
    if api_key:
        headers["Authorization"] = "Bearer " + api_key
    request = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            payload = json.load(response)
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", "replace")[:500]
        raise RuntimeError(f"OpenAlex returned HTTP {exc.code}: {detail}") from exc
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"could not read OpenAlex: {exc}") from exc
    results = payload.get("results")
    if not isinstance(results, list):
        raise RuntimeError("OpenAlex response has no results list")
    return results


def fetch(dois: list[str], api_key: str = "") -> dict[str, tuple[int, str]]:
    out: dict[str, tuple[int, str]] = {}
    for group in batches(dois):
        for work in request_batch(group, api_key):
            doi_url = work.get("doi") or ""
            doi = doi_url.removeprefix("https://doi.org/").lower()
            count, work_id = work.get("cited_by_count"), work.get("id") or ""
            if doi not in group or not isinstance(count, int) or count < 0 or not work_id:
                raise RuntimeError(f"unexpected OpenAlex work record: {work!r}")
            if doi in out:
                raise RuntimeError(f"OpenAlex returned DOI {doi!r} more than once")
            out[doi] = (count, work_id)
    return out


def write_snapshot(path: Path, dois: list[str], found: dict[str, tuple[int, str]], day: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", encoding="utf-8", newline="", dir=path.parent,
                                     prefix=path.name + ".", delete=False) as fh:
        tmp = Path(fh.name)
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for doi in dois:
            count, work_id = found.get(doi, ("", ""))
            writer.writerow({"doi": doi, "citations": count, "openalex_id": work_id,
                             "retrieved": day})
    tmp.replace(path)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=TSV,
                    help=f"snapshot path (default: {TSV.relative_to(ROOT)})")
    args = ap.parse_args(argv)
    dois = sorted({r["doi"].strip().lower() for r in catalogue.ours() if r["doi"].strip()})
    try:
        found = fetch(dois, os.environ.get("OPENALEX_API_KEY", ""))
    except RuntimeError as exc:
        print(f"fetch_citations: {exc}", file=sys.stderr)
        return 1
    write_snapshot(args.out, dois, found, date.today().isoformat())
    missing = [d for d in dois if d not in found]
    print(f"wrote {args.out}: {len(found)} citation counts, {len(missing)} unresolved DOIs")
    for doi in missing:
        print(f"  unresolved  {doi}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
