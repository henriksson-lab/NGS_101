"""
Read access to the protocol catalogue, for the notes and for future tooling.

The table itself is built by fetch_scg_lib_structs.py and is the only copy of these facts.
Anything that wants to quote a count, a protocol list or a directory name reads it from
here -- including the Markdown notes, through `lib/mdfacts.py`:

    {{= catalogue.n_protocols() }} protocols, {{= catalogue.n_papers() }} distinct papers

One row per (protocol, paper) pair, because the relation is many-to-many in both
directions: SMART-seq family cites six papers, and Smart-seq2's paper is the RNA half of
scM&T-seq, scMT-seq and scNMT-seq as well. `role` separates a method's own paper(s) from
the incidental technique references its page also links.
"""

from __future__ import annotations

import csv
from functools import lru_cache
from pathlib import Path

TSV = Path(__file__).resolve().parents[1] / "scg_lib_structs.tsv"
OURS = Path(__file__).resolve().parents[1] / "ours.tsv"
SOURCE_REPO = "https://github.com/Teichlab/scg_lib_structs"
SOURCE_PAGES = "https://teichlab.github.io/scg_lib_structs/"

# ours.tsv's columns, and the closed vocabularies of three of them.
OURS_COLUMNS = ["dir", "protocol", "doi", "status", "section", "modality", "note"]
STATUSES = ("documented", "draft", "notes")
# Which part of the website a directory is listed in: `published` protocols go into the
# searchable list, `wip` (our own unpublished designs) into a separate section.
SECTIONS = ("published", "wip")
MODALITIES = ("DNA", "RNA", "multi")
# What each scg_lib_structs category implies about the molecule sequenced. The other two
# categories ("TODO list", "Ours, not in scg_lib_structs") imply nothing, which is why
# `modality` is recorded in ours.tsv rather than derived.
CATEGORY_MODALITY = {
    "Gene expression": "RNA",
    "Chromatin accessibility and protein-DNA interactions": "DNA",
    "Genomic DNA or DNA methylation": "DNA",
    "Multi-Omics": "multi",
}


@lru_cache(maxsize=1)
def rows() -> tuple[dict, ...]:
    with TSV.open(encoding="utf-8") as fh:
        return tuple(csv.DictReader(fh, delimiter="\t"))


@lru_cache(maxsize=1)
def ours() -> tuple[dict, ...]:
    """What THIS repo has tackled: one row per protocol directory (catalogue/ours.tsv).

    This is the only place our own coverage is recorded. The directory name is derived
    from the same scheme the catalogue uses, so `dir` is both the directory on disk and
    the join key back into the scraped table (by DOI).
    """
    with OURS.open(encoding="utf-8") as fh:
        return tuple(csv.DictReader(fh, delimiter="\t"))


def ours_in(section: str) -> list[dict]:
    """Our directories listed in one website section (`published` or `wip`)."""
    if section not in SECTIONS:
        raise ValueError(f"unknown section {section!r}; one of {SECTIONS}")
    return [r for r in ours() if r["section"] == section]


def published() -> list[dict]:
    return ours_in("published")


def wip() -> list[dict]:
    return ours_in("wip")


def rows_for_dir(d: str) -> list[dict]:
    """Every catalogue row joined to one of our directories (all its papers)."""
    return [r for r in rows() if r["our_dir"] == d]


def ours_by_doi() -> dict[str, dict]:
    return {r["doi"]: r for r in ours() if r["doi"]}


def covered() -> list[dict]:
    """Rows of the table that this repo documents (its `our_dir` is filled in)."""
    return [r for r in rows() if r["our_dir"]]


def overlap() -> list[dict]:
    """Rows we cover that upstream also lists -- i.e. where we can cross-check ourselves."""
    return [r for r in rows() if r["our_dir"] and r["source"] == "scg_lib_structs"]


def untackled() -> list[dict]:
    """Defining rows for protocols nobody here has documented yet: the worklist."""
    return [r for r in rows() if r["is_defining"] == "yes" and not r["our_dir"]]


def n_ours() -> int:
    return len(ours())


def n_ours_documented() -> int:
    return sum(1 for r in ours() if r["status"] == "documented")


def ours_table() -> str:
    """A Markdown table of our own directories -- generated, never hand-listed."""
    head = "| protocol | directory | paper | status | section |"
    rule = "|---|---|---|---|---|"
    body = []
    for r in ours():
        doi = f"[{r['doi']}](https://doi.org/{r['doi']})" if r["doi"] else "&mdash;"
        body.append(f"| {r['protocol']} | `{r['dir']}/` | {doi} | {r['status']} "
                    f"| {r['section']} |")
    return "\n".join([head, rule, *body])


def documented() -> list[dict]:
    """Rows for protocols scg_lib_structs has already drawn a page for."""
    return [r for r in rows() if r["documented"] == "yes"]


def from_scg() -> list[dict]:
    return [r for r in rows() if r["source"] == "scg_lib_structs"]


def todo() -> list[dict]:
    """Rows from its TODO list: a method named, with a paper, but no page yet."""
    return [r for r in rows() if r["documented"] == "todo"]


def primary() -> list[dict]:
    """Rows that are a method's own paper, not an incidental citation."""
    return [r for r in rows() if r["role"] == "primary"]


def defining() -> list[dict]:
    """One row per protocol: the paper its directory is named after."""
    return [r for r in rows() if r["is_defining"] == "yes"]


def protocols() -> list[str]:
    return sorted({r["protocol"] for r in rows()})


def categories() -> list[str]:
    seen: list[str] = []
    for r in rows():
        if r["category"] not in seen:
            seen.append(r["category"])
    return seen


def n_protocols() -> int:
    return len(protocols())


def n_documented() -> int:
    return len({r["protocol"] for r in documented()})


def n_todo() -> int:
    return len({r["protocol"] for r in todo()})


def n_papers() -> int:
    """Distinct papers, counting a DOI once however many protocols cite it."""
    return len({r["doi"] for r in rows() if r["doi"]})


def n_rows() -> int:
    return len(rows())


def with_doi() -> int:
    return sum(1 for r in rows() if r["doi"])


def with_pmid() -> int:
    return sum(1 for r in rows() if r["pmid"])


def by_category() -> list[tuple[str, int, int]]:
    """[(category, n protocols, n rows)] in README order."""
    out = []
    for c in categories():
        rs = [r for r in rows() if r["category"] == c]
        out.append((c, len({r["protocol"] for r in rs}), len(rs)))
    return out


def multi_paper() -> list[tuple[str, int]]:
    """Protocols whose own papers number more than one, commonest first."""
    counts: dict[str, int] = {}
    for r in primary():
        counts[r["protocol"]] = counts.get(r["protocol"], 0) + 1
    return sorted(((p, n) for p, n in counts.items() if n > 1),
                  key=lambda x: (-x[1], x[0]))


def shared_papers() -> list[tuple[str, str, list[str]]]:
    """[(doi, title, [protocols])] for papers that are the own-paper of several methods."""
    by: dict[str, set[str]] = {}
    title: dict[str, str] = {}
    for r in primary():
        if not r["doi"]:
            continue
        by.setdefault(r["doi"], set()).add(r["protocol"])
        title[r["doi"]] = r["title"]
    return sorted(((d, title[d], sorted(p)) for d, p in by.items() if len(p) > 1),
                  key=lambda x: (-len(x[2]), x[0]))


def year_range() -> tuple[str, str]:
    ys = sorted(r["year"] for r in rows() if r["year"])
    return (ys[0], ys[-1]) if ys else ("", "")


def markdown_table(rs: list[dict], columns=("protocol", "year", "journal", "doi")) -> str:
    """A Markdown table of `rs` -- so a note never hand-types catalogue contents."""
    head = "| " + " | ".join(columns) + " |"
    rule = "|" + "|".join("---" for _ in columns) + "|"
    body = []
    for r in rs:
        cells = []
        for c in columns:
            v = r.get(c, "")
            if c == "doi" and v:
                v = f"[{v}](https://doi.org/{v})"
            elif c == "protocol" and r.get("scg_page"):
                v = f"[{v}]({r['scg_page']})"
            cells.append(v)
        body.append("| " + " | ".join(cells) + " |")
    return "\n".join([head, rule, *body])
