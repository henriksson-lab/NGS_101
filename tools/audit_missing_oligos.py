#!/usr/bin/env python3
"""Inventory rendered protocol pages that declare unresolved oligo sequence.

This is a discovery and review ledger, not a sequence oracle. It extracts explicit
limitations from the same schematic pages users see, then preserves a curator's source
review beside the current signal digest. If page wording or coverage changes after a
review, the digest changes and the row becomes visibly stale.

Run after building protocol pages::

    python3 tools/audit_missing_oligos.py
    python3 tools/audit_missing_oligos.py --write

The second form updates ``catalogue/oligo_source_audit.tsv`` while retaining completed
review fields. Third-party documents are never copied into the repository.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import re
import sys
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "catalogue" / "tools"))

import build_index  # noqa: E402
import catalogue as cat  # noqa: E402

OUT = ROOT / "catalogue" / "oligo_source_audit.tsv"
COLUMNS = ("dir", "protocol", "doi", "active", "signal_digest", "page_signals",
           "audit_status", "sources_checked", "conclusion", "reviewed_on",
           "reviewed_digest")
STATUSES = ("pending", "resolved", "confirmed_missing", "source_unavailable",
            "commercial_proprietary", "mixed", "not_applicable")

# These phrases claim uncertainty or absence. Requiring an oligo-related object keeps
# generic unknown biological outcomes out of the inventory.
LIMITATION = re.compile(
    r"(?i)\b(?:unavailable|proprietary|undisclosed|unpublished|unknown|unresolved|"
    r"omitted|missing|not (?:public|published|disclosed|available|given|printed|"
    r"reported|specified|provided|known|represented|resolved|verified|transcribed)|"
    r"does not (?:publish|print|disclose|expose|give|report|specify|provide)|"
    r"cannot be (?:drawn|computed|resolved|verified)|role[- ]only|shown by role|"
    r"source limitation|sequence limitation)\b")
OLIGO_OBJECT = re.compile(
    r"(?i)\b(?:primer|adapter|adaptor|linker|oligo|barcode|sequence|overhang|handle|"
    r"transposome|whitelist)s?\b")


class LimitationBlocks(HTMLParser):
    """Collect human-readable blocks, excluding scripts and the graph audit dump."""

    def __init__(self) -> None:
        super().__init__()
        self.stack: list[tuple[str, list[str]]] = []
        self.blocks: list[str] = []
        self.skip = 0

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attributes = {k: v or "" for k, v in attrs}
        classes = set(attributes.get("class", "").split())
        if tag in ("script", "style"):
            self.skip += 1
        wanted = tag in ("figcaption", "small") or (
            tag == "p" and "protocol-completeness" not in classes) or (
            tag == "div" and classes & {"caveat", "warning", "notice"})
        if wanted:
            self.stack.append((tag, []))
        if tag == "br" and not self.skip:
            for _, chunks in self.stack:
                chunks.append(" ")

    def handle_data(self, data: str) -> None:
        if not self.skip:
            for _, chunks in self.stack:
                chunks.append(data)

    def handle_endtag(self, tag: str) -> None:
        if tag in ("script", "style") and self.skip:
            self.skip -= 1
        for i in range(len(self.stack) - 1, -1, -1):
            if self.stack[i][0] == tag:
                _, chunks = self.stack.pop(i)
                text = " ".join("".join(chunks).split())
                if text:
                    self.blocks.append(text)
                break


def limitations(page: Path) -> tuple[str, ...]:
    parser = LimitationBlocks()
    parser.feed(page.read_text(encoding="utf-8", errors="replace"))
    found: list[str] = []
    for block in parser.blocks:
        # A whole warning box often contrasts an exact oligo statement with an
        # unrelated unavailable Methods detail.  Match within the same prose clause,
        # otherwise those pages become false source-review obligations.
        clauses = re.split(r"(?<=[.!?;])\s+|\s+[—–]\s+", block)
        for clause in clauses:
            clause = clause.strip()
            # Workflow captions name discarded branches such as "unresolved GEM
            # products".  That describes sample disposition, not an unresolved
            # oligo; nested figure text can otherwise place a later barcode word in
            # the same HTML-parser block and create a false primer-audit candidate.
            if clause.lower().startswith("other outputs:"):
                continue
            if (LIMITATION.search(clause) and OLIGO_OBJECT.search(clause)
                    and clause not in found):
                found.append(clause)
    return tuple(found)


def digest(signals: tuple[str, ...]) -> str:
    canonical = "\n".join(signals).encode("utf-8")
    return hashlib.sha256(canonical).hexdigest()[:16]


def candidates() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for record in cat.ours():
        if record["status"] == "notes":
            continue
        page = build_index.diagram_out(record["dir"])
        if page is None or not page.exists():
            continue
        signals = limitations(page)
        if not signals:
            continue
        rows.append({
            "dir": record["dir"],
            "protocol": record["protocol"],
            "doi": record["doi"],
            "active": "yes",
            "signal_digest": digest(signals),
            "page_signals": " || ".join(s.replace("\t", " ") for s in signals),
            "audit_status": "pending",
            "sources_checked": "",
            "conclusion": "",
            "reviewed_on": "",
            "reviewed_digest": "",
        })
    return rows


def existing() -> dict[str, dict[str, str]]:
    if not OUT.exists():
        return {}
    with OUT.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    rows = [{field: (row.get(field) or "") for field in COLUMNS} for row in rows]
    if rows and "active" not in rows[0]:
        for row in rows:
            row["active"] = "yes"
    if rows and set(rows[0]) != set(COLUMNS):
        raise ValueError(f"{OUT}: columns are {tuple(rows[0])!r}, expected {COLUMNS!r}")
    return {row["dir"]: row for row in rows}


def merged() -> list[dict[str, str]]:
    old = existing()
    rows = candidates()
    current = {row["dir"] for row in rows}
    for row in rows:
        before = old.get(row["dir"])
        if before:
            for field in ("audit_status", "sources_checked", "conclusion",
                          "reviewed_on", "reviewed_digest"):
                row[field] = before[field]
        if row["audit_status"] not in STATUSES:
            raise ValueError(f"{row['dir']}: unknown audit status {row['audit_status']!r}")
    for directory, before in old.items():
        # An unreviewed discovery row has no historical evidentiary value once a
        # more precise extractor proves it was a false positive.  Completed reviews
        # remain visible even after the page no longer declares the gap.
        if directory not in current and before["audit_status"] != "pending":
            before["active"] = "no"
            rows.append(before)
    order = {row["dir"]: i for i, row in enumerate(cat.ours())}
    rows.sort(key=lambda row: order.get(row["dir"], len(order)))
    return rows


def write(rows: list[dict[str, str]]) -> None:
    with OUT.open("w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            values = [row.get(field, "").replace("\t", " ").replace("\n", " ")
                      for field in COLUMNS]
            # Pending rows deliberately have empty review columns. Omitting their
            # trailing separators keeps the checked-in TSV clean; DictReader restores
            # those absent trailing cells to empty strings on the next read.
            handle.write("\t".join(values).rstrip("\t") + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--write", action="store_true", help=f"update {OUT.relative_to(ROOT)}")
    parser.add_argument("--review", metavar="DIR", help="record a completed review")
    parser.add_argument("--status", choices=STATUSES)
    parser.add_argument("--sources", help="authoritative materials actually read")
    parser.add_argument("--conclusion", help="what was recovered or remains unavailable")
    args = parser.parse_args()
    rows = merged()
    if args.review:
        if not (args.status and args.sources and args.conclusion):
            parser.error("--review requires --status, --sources and --conclusion")
        matches = [r for r in rows if r["dir"] == args.review]
        if not matches and args.status in ("resolved", "not_applicable"):
            records = [r for r in cat.ours() if r["dir"] == args.review]
            if len(records) == 1:
                record = records[0]
                row = {field: "" for field in COLUMNS}
                row.update(dir=record["dir"], protocol=record["protocol"],
                           doi=record["doi"], active="no")
                rows.append(row)
                matches = [row]
        if len(matches) != 1:
            parser.error(f"{args.review!r} is not exactly one inventory row")
        row = matches[0]
        row["audit_status"] = args.status
        row["sources_checked"] = args.sources.replace("\t", " ")
        row["conclusion"] = args.conclusion.replace("\t", " ")
        row["reviewed_on"] = dt.date.today().isoformat()
        row["reviewed_digest"] = row["signal_digest"]
        write(rows)
    elif args.write:
        write(rows)

    stale = [r for r in rows if r["active"] == "yes" and r["audit_status"] != "pending"
             and r["reviewed_digest"] != r["signal_digest"]]
    pending = [r for r in rows if r["active"] == "yes" and r["audit_status"] == "pending"]
    active = sum(r["active"] == "yes" for r in rows)
    print(f"{active} protocol pages declare a possible unresolved oligo; "
          f"{len(pending)} pending source review, {len(stale)} reviewed rows stale, "
          "0 cleared signals awaiting disposition")
    for row in stale:
        print(f"STALE\t{row['dir']}\t{row['signal_digest']}\t{row['reviewed_digest']}")
    return 1 if stale else 0


if __name__ == "__main__":
    raise SystemExit(main())
