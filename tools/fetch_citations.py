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
        # OpenAlex accepts bare DOI values in an OR filter. Supplying full doi.org URLs
        # makes the embedded ':' ambiguous to its filter parser and can return unrelated
        # works, so keep the transport URL out of the identifier values.
        "filter": "doi:" + "|".join(dois),
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


def request_one(doi: str, api_key: str = "") -> dict | None:
    """Resolve one requested DOI, preserving that DOI even when OpenAlex canonicalizes it.

    A preprint DOI can resolve to a work whose returned ``doi`` is the journal DOI. Batch
    responses then cannot say which requested alias produced the record, so unresolved
    batch inputs are retried through the identifier endpoint one by one.
    """
    ident = urllib.parse.quote("https://doi.org/" + doi, safe="")
    url = API + "/" + ident + "?select=id,doi,cited_by_count"
    headers = {"User-Agent": USER_AGENT}
    if api_key:
        headers["Authorization"] = "Bearer " + api_key
    try:
        with urllib.request.urlopen(urllib.request.Request(url, headers=headers),
                                    timeout=60) as response:
            return json.load(response)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return None
        detail = exc.read().decode("utf-8", "replace")[:500]
        raise RuntimeError(f"OpenAlex returned HTTP {exc.code}: {detail}") from exc
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"could not read OpenAlex: {exc}") from exc


def fetch(dois: list[str], api_key: str = "") -> dict[str, tuple[int, str]]:
    out: dict[str, tuple[int, str]] = {}
    ambiguous: set[str] = set()
    for group in batches(dois):
        for work in request_batch(group, api_key):
            doi_url = work.get("doi") or ""
            doi = doi_url.removeprefix("https://doi.org/").lower()
            count, work_id = work.get("cited_by_count"), work.get("id") or ""
            if doi not in group:
                continue                    # canonical record returned for a DOI alias
            if not isinstance(count, int) or count < 0 or not work_id:
                raise RuntimeError(f"unexpected OpenAlex work record: {work!r}")
            candidate = (count, work_id)
            if doi in ambiguous:
                continue
            if doi in out and out[doi] != candidate:
                # OpenAlex can retain multiple works for one DOI (e.g. nbt.2137).
                # Resolve through its DOI endpoint below, never by batch order or
                # by choosing whichever record happens to have more citations.
                del out[doi]
                ambiguous.add(doi)
            else:
                out[doi] = candidate
    for doi in (d for d in dois if d not in out):
        work = request_one(doi, api_key)
        if work is None:
            continue
        count, work_id = work.get("cited_by_count"), work.get("id") or ""
        if not isinstance(count, int) or count < 0 or not work_id:
            raise RuntimeError(f"unexpected OpenAlex work record: {work!r}")
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
