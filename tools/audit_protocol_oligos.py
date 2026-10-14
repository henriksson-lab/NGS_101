#!/usr/bin/env python3
"""Track source-level oligo review for every protocol in catalogue/ours.tsv.

The rendered-gap ledger in audit_missing_oligos.py finds explicit uncertainty claims.
This separate census also covers pages that currently claim no gap and notes-only
protocols. A review is evidence of a source check, not proof of sequence completeness.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "catalogue" / "tools"))
import catalogue as cat  # noqa: E402

OUT = ROOT / "catalogue" / "protocol_oligo_audit.tsv"
COLUMNS = ("dir", "protocol", "doi", "audit_status", "sources_checked",
           "oligo_roles_checked", "open_roles", "conclusion", "reviewed_on")
STATUSES = ("pending", "reviewed", "limited", "not_applicable")


def rows() -> list[dict[str, str]]:
    previous: dict[str, dict[str, str]] = {}
    if OUT.exists():
        with OUT.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != COLUMNS:
                raise ValueError(f"{OUT}: unexpected columns {reader.fieldnames!r}")
            previous = {row["dir"]: row for row in reader}
    result = []
    for record in cat.ours():
        row = {field: "" for field in COLUMNS}
        row.update({field: record[field] for field in ("dir", "protocol", "doi")})
        row["audit_status"] = "pending"
        if record["dir"] in previous:
            for field in COLUMNS[3:]:
                row[field] = previous[record["dir"]].get(field) or ""
        if row["audit_status"] not in STATUSES:
            raise ValueError(f"{row['dir']}: unknown status {row['audit_status']!r}")
        result.append(row)
    return result


def write(items: list[dict[str, str]]) -> None:
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for row in items:
            writer.writerow({key: row[key].replace("\t", " ").replace("\n", " ")
                             for key in COLUMNS})


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--write", action="store_true", help="sync rows to catalogue")
    parser.add_argument("--review", metavar="DIR", help="record one protocol review")
    parser.add_argument("--status", choices=STATUSES)
    parser.add_argument("--sources", help="authoritative sources actually inspected")
    parser.add_argument("--roles", help="oligo roles checked against those sources")
    parser.add_argument("--open-roles", default="", help="roles unresolved after review")
    parser.add_argument("--conclusion", help="short evidence boundary and any inference")
    args = parser.parse_args()
    items = rows()
    if args.review:
        if args.status not in STATUSES[1:] or not all((args.sources, args.roles, args.conclusion)):
            parser.error("--review requires a non-pending --status, --sources, --roles and --conclusion")
        match = next((row for row in items if row["dir"] == args.review), None)
        if match is None:
            parser.error(f"unknown protocol directory: {args.review}")
        match.update(audit_status=args.status, sources_checked=args.sources,
                     oligo_roles_checked=args.roles, open_roles=args.open_roles,
                     conclusion=args.conclusion, reviewed_on=dt.date.today().isoformat())
        write(items)
    elif args.write:
        write(items)
    counts = Counter(row["audit_status"] for row in items)
    print(f"{len(items)} protocols: " + ", ".join(f"{s} {counts[s]}" for s in STATUSES))
    if not OUT.exists() or len(items) != sum(1 for _ in OUT.open(encoding="utf-8")) - 1:
        print("Run with --write to sync the catalogue census.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
