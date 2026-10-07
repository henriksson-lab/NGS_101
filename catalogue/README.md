# The protocol catalogue

A worklist for documenting NGS library chemistries base-by-base, in the style the rest of
this repo uses. The aim is near-complete coverage, so the first job is knowing what exists.

**Nothing here is typed by hand.** `tools/fetch_scg_lib_structs.py` scrapes
[Teichlab/scg_lib_structs](https://github.com/Teichlab/scg_lib_structs) — the reference
collection this repo's house style is modelled on — and resolves every paper it finds
against Crossref and NCBI. Re-running it updates the table; the numbers in this note come
from the table itself (via `lib/mdfacts.py`), so they cannot drift from it.

## What is in it

`scg_lib_structs.tsv`: **{{= catalogue.n_rows() }} rows**
({{= len(catalogue.from_scg()) }} scraped from upstream, {{= catalogue.n_rows() - len(catalogue.from_scg()) }} added for
protocols only we cover), covering
**{{= catalogue.n_protocols() }} protocols** and **{{= catalogue.n_papers() }} distinct
papers**, published {{= catalogue.year_range()[0] }}–{{= catalogue.year_range()[1] }}.

- **{{= catalogue.n_documented() }} protocols** already have a drawn page upstream
  (`documented = yes`), which is both a cross-check and a style reference for us.
- **{{= catalogue.n_todo() }} protocols** are on the upstream TODO list
  (`documented = todo`): named, with a paper, but never drawn. These are open ground.

| category | protocols | rows |
|---|---|---|
{{= "\n".join(f"| {c} | {np} | {nr} |" for c, np, nr in catalogue.by_category()) }}

Coverage of identifiers: **{{= catalogue.with_doi() }}/{{= catalogue.n_rows() }} rows have a
DOI**, {{= catalogue.with_pmid() }} have a PMID. The gaps are real and marked in `note`:
four vendor protocols documented from a kit manual with no publication at all, and one
Protocol Exchange entry whose DOI no metadata service confirms — recorded as absent rather
than guessed.

## One row per (protocol, paper), because it is not one-to-one

Both directions are many-to-many, which is why `protocol` is not a key:

- **A protocol described in several papers.** {{= len(catalogue.multi_paper()) }} of them,
  led by {{= "; ".join(f"{p} ({n})" for p, n in catalogue.multi_paper()[:3]) }}.
- **A paper that defines several protocols.** The real cases are method compositions:

{{= catalogue.markdown_table([{"protocol": ", ".join(ps), "year": "", "journal": t[:60], "doi": d} for d, t, ps in catalogue.shared_papers()[:5]], ("protocol", "journal", "doi")) }}

`role` separates a method's own paper from the incidental technique references its page
also links — without that split, the 2010 semi-suppressive-PCR paper looks like the source
of 20 unrelated protocols, because 20 pages cite it in passing. `is_defining` marks the one
paper per protocol that names its directory.

## Directory naming

```
<protocol-slug>__<doi, with every "/" replaced by "+">

smart-seq-family__10.1038+nbt.2282
sci-rna-seq-family__10.1126+science.aam8940
limca__10.21203+rs.3.rs-3210240+v1
```

- **DOI, not PMID.** Every preprint has a DOI the day it is posted, and in this table
  {{= catalogue.with_doi() }} rows carry one against {{= catalogue.with_pmid() }} with a
  PMID. The PMID is kept as a column for anyone who wants it.
- **Reversible.** `+` is the only substitution, it is filesystem- and shell-safe, and no
  DOI contains one (the selftest checks every DOI in the table). Replacing only the first
  `/` would not do — a Research Square DOI has two — and `_` would be ambiguous, because
  DOI suffixes contain underscores themselves (`10.1038/nbt0400_424`). `__` separates the
  name from the DOI; a slugified name never contains a run of underscores.
- **One directory per protocol, named for its *defining* paper** — the earliest paper its
  page cites as its own. A method described in several papers still gets one directory; the
  other papers stay in the table as further rows.
- **A protocol with no paper** (a vendor kit) is named by its name alone.

## Using it

```sh
python3 catalogue/tools/fetch_scg_lib_structs.py     # rebuild the table (cached, re-runnable)
python3 catalogue/tools/fetch_scg_lib_structs.py --offline   # from the cache only
python3 catalogue/tools/selftest.py                  # check the table's invariants
```

The download cache holds third-party HTML. It is read as text, parsed with regexes, never
executed, and kept outside the repo.

## Status

This is the worklist, not a plan of record. Still to decide: which protocols are worth a
full base-by-base page, in what order, and what a second source (beyond scg_lib_structs)
should be — the TODO list is upstream's backlog, not necessarily ours, and the catalogue
covers single-cell genomics only. Bulk protocols, long-read chemistries and spatial
methods are not in it yet.

## What we have tackled

Recorded in `ours.tsv`, which is also the source of the directory names on disk and of the
`our_dir` column in the scraped table — so the two cannot disagree.
**{{= catalogue.n_ours_documented() }} documented, {{= catalogue.n_ours() }} in total:**

{{= catalogue.ours_table() }}

The join is by **DOI**, not by name: a protocol's name differs between upstream's list and
ours, but its defining paper does not. {{= len(catalogue.overlap()) }} rows overlap with
upstream (where we can check ourselves against their drawing); the rest are protocols
upstream does not cover — small-RNA, pooled-screen readouts, and our own unpublished
design — and are appended to the table with `source = ours`, so it is the whole worklist
rather than only upstream's part of it.

**{{= len(catalogue.untackled()) }} protocols are still untouched here.**
