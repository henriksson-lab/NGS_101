#!/usr/bin/env python3
"""Report and validate independent source reviews for every catalogue protocol.

This is an acceptance gate, not a source-reading substitute. A historical
``reviewed`` row in either oligo ledger does not satisfy the two-pass rule.
Reviews live in catalogue/source_reviews/<protocol-dir>.json; see its README.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REVIEWS = ROOT / "catalogue" / "source_reviews"
CATALOGUE = ROOT / "catalogue" / "ours.tsv"
STATUSES = {"in_progress", "blocked", "reviewed", "accepted_missing"}
ROLES = {"defining", "supplement", "detailed_protocol", "vendor",
         "predecessor", "linked_material", "other"}


def catalogue() -> dict[str, dict[str, str]]:
    with CATALOGUE.open(encoding="utf-8", newline="") as handle:
        return {row["dir"]: row for row in csv.DictReader(handle, delimiter="\t")}


def _filled(obj: dict, *keys: str) -> bool:
    return all(isinstance(obj.get(key), str) and obj[key].strip() for key in keys)


def _pass_ok(value: object) -> bool:
    return (isinstance(value, dict)
            and _filled(value, "reviewer", "date", "locations", "finding"))


def _two_passes(obj: dict, label: str, issues: list[str]) -> None:
    first, second = obj.get("pass_1"), obj.get("pass_2")
    if not _pass_ok(first) or not _pass_ok(second):
        issues.append(label + ": two located source checks are required")
    elif first["reviewer"].strip().casefold() == second["reviewer"].strip().casefold():
        issues.append(label + ": second check must have a different reviewer")


def validate(name: str, review: object, known: dict[str, dict[str, str]]) -> list[str]:
    issues: list[str] = []
    if not isinstance(review, dict):
        return [name + ": review must be a JSON object"]
    if name not in known or review.get("dir") != name:
        issues.append(name + ": directory is absent from catalogue or disagrees with file")
    status = review.get("status")
    if status not in STATUSES:
        issues.append(name + ": invalid review status")
        return issues
    sources = review.get("source_checklist")
    if not isinstance(sources, list) or not sources:
        issues.append(name + ": source checklist is empty")
        sources = []
    ids: set[str] = set()
    for i, source in enumerate(sources, 1):
        label = f"{name}: source {i}"
        if not isinstance(source, dict):
            issues.append(label + " is not an object")
            continue
        if not _filled(source, "id", "role", "locator", "availability", "relevance"):
            issues.append(label + ": id, role, locator, availability, relevance required")
        if source.get("id") in ids:
            issues.append(label + ": duplicate source ID")
        ids.add(source.get("id", ""))
        if source.get("role") not in ROLES:
            issues.append(label + ": unknown source role")
        if source.get("availability") not in {"read", "unavailable"}:
            issues.append(label + ": availability must be read or unavailable")
        if source.get("relevance") not in {"relevant", "not_relevant"}:
            issues.append(label + ": relevance must be relevant or not_relevant")
        if source.get("role") == "linked_material" and not source.get("discovered_from"):
            issues.append(label + ": linked material needs the referring source ID")
        if source.get("availability") == "unavailable" and not source.get("reason"):
            issues.append(label + ": unavailable source needs a reason")
        if source.get("relevance") == "not_relevant" and not source.get("reason"):
            issues.append(label + ": excluded source needs a reason")
        if status in {"reviewed", "accepted_missing"}:
            if source.get("relevance") == "relevant" and source.get("availability") != "read":
                issues.append(label + ": relevant source unavailable; review is blocked")
            _two_passes(source, label, issues)
    if sources and not any(s.get("role") == "defining" for s in sources if isinstance(s, dict)):
        issues.append(name + ": checklist lacks a defining article or vendor source")
    links = review.get("linked_information")
    if not isinstance(links, dict) or links.get("result") not in {
            "enumerated", "none_found", "blocked"}:
        issues.append(name + ": linked_information needs a result")
    elif status in {"reviewed", "accepted_missing"}:
        if links["result"] == "blocked":
            issues.append(name + ": linked information remains blocked")
        _two_passes(links, name + ": linked information", issues)
        if links["result"] == "enumerated" and not any(
                s.get("role") == "linked_material" for s in sources if isinstance(s, dict)):
            issues.append(name + ": enumerated links need checklist entries")
    claims = review.get("missing_claims")
    if not isinstance(claims, list):
        issues.append(name + ": missing_claims must be a list")
        claims = []
    for i, claim in enumerate(claims, 1):
        if not isinstance(claim, dict) or not _filled(
                claim, "claim", "disposition", "evidence"):
            issues.append(f"{name}: missing claim {i} lacks claim, disposition or evidence")
        elif claim["disposition"] not in {"resolved", "confirmed_missing", "blocked"}:
            issues.append(f"{name}: missing claim {i} has unknown disposition")
    if status in {"reviewed", "accepted_missing"} and any(
            isinstance(c, dict) and c.get("disposition") == "blocked" for c in claims):
        issues.append(name + ": unresolved missing claim blocks completed review")
    if status == "reviewed" and any(
            isinstance(c, dict) and c.get("disposition") == "confirmed_missing" for c in claims):
        issues.append(name + ": confirmed missing claim needs accepted_missing status")
    if status == "accepted_missing":
        if not any(isinstance(c, dict) and c.get("disposition") == "confirmed_missing"
                   for c in claims):
            issues.append(name + ": accepted_missing needs a confirmed missing claim")
        if not _filled(review, "conclusion", "public_page", "public_checklist_marker"):
            issues.append(name + ": acceptance needs conclusion and public checklist marker")
        else:
            page = ROOT / review["public_page"]
            if not page.is_file() or review["public_checklist_marker"] not in page.read_text(
                    encoding="utf-8", errors="replace"):
                issues.append(name + ": public page does not display the source checklist")
    return issues


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", metavar="DIR", help="validate one protocol and fail if absent")
    parser.add_argument("--pending", action="store_true", help="list protocols lacking a review")
    args = parser.parse_args()
    known = catalogue()
    if args.check and args.check not in known:
        parser.error("unknown protocol directory: " + args.check)
    names = [args.check] if args.check else list(known)
    counts: Counter[str] = Counter()
    failures: list[str] = []
    for name in names:
        path = REVIEWS / (name + ".json")
        if not path.is_file():
            counts["pending"] += 1
            if args.check:
                failures.append(name + ": no source review file")
            elif args.pending:
                print(name)
            continue
        try:
            review = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as exc:
            failures.append(name + ": " + str(exc))
            continue
        counts[review.get("status", "invalid")] += 1
        failures.extend(validate(name, review, known))
    extra = {p.stem for p in REVIEWS.glob("*.json")} - set(known)
    failures.extend(name + ": review has no catalogue protocol" for name in sorted(extra))
    print(f"{len(known)} protocols: " + ", ".join(f"{key} {counts[key]}" for key in
          ("accepted_missing", "reviewed", "in_progress", "blocked", "pending")))
    for problem in failures:
        print("FAIL", problem)
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
