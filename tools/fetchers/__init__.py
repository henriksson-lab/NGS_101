"""Fetchers: special-case source handlers that tools/get_sources.py consults.

get_sources.py knows the generic routes for a paper (PMC / the PMC Cloud bucket, Europe
PMC, bioRxiv, Springer supplements, preprint <-> published links). Everything that needs
special handling -- a vendor's CDN, a GitHub repository, protocols.io, eLife's API, the
next publisher -- is a *fetcher*: a module in this package that registers a handler.

Two ways a handler gets work:

1. **Automatically**, by pattern. A handler can claim
     doi=r"^10\\.17504/protocols\\.io\\."   any paper DOI of the protocol matching this
     url=r"^https://www\\.protocols\\.io/"  any extra-source URL matching this
     auto=lambda job: ...                  a predicate on the protocol (e.g. "has an
                                           scg_lib_structs page")
2. **Declaratively**, from catalogue/extra_sources.tsv -- one row per extra source a
   protocol needs (tab-separated; '#' starts a comment line):

     slug                         source                                    what   why
     10x-chromium-single-cell-5-vdj  https://cdn.10xgenomics.com/...CG000331...pdf  10x user guide CG000331 Rev E  primary oligo source; no DOI
     bd-rhapsody__10.1038+...     scg:data/BD/*                             BD user guides  vendor kit
     scrb-seq-mcscrb-seq__...     protocolsio:10.17504/protocols.io.nrkdd4w mcSCRB-seq protocol  unblocked TSO
     strt-seq-family__...         doi:10.1038/s41598-017-16546-4            STRT-seq-2i journal version  only the preprint is catalogued

   `source` is either a URL (routed to the handler whose `url` pattern matches, else the
   generic `url` handler, which downloads and validates it) or `NAME:ARG` for a
   registered handler NAME (`doi:` and `pmc:` are built in to get_sources: fetch that
   paper the normal way).

Adding a fetcher
----------------
Create tools/fetchers/<name>.py:

    from . import register
    from .base import say, get_json

    @register("myname", doi=r"^10\\.1234/", url=r"^https://example\\.org/",
              doc="one line: what it fetches and from where")
    def fetch(job, arg: str, what: str = "") -> None:
        # job.folder.save(filename, url, what) downloads, validates, converts, records;
        # a failed/invalid download becomes a '(manual)' manifest row by itself.
        # job.paper(doi) fetches another paper the generic way; job.rows are the
        # protocol's catalogue rows; job.slug its directory name.
        ...

then add its module name to MODULES below. Handlers must be idempotent (Folder.save skips
what is on disk), must not write outside job.folder.path, and should record a
`(manual)` row (job.folder.manual) for anything they know exists but cannot fetch. Add a
row to catalogue/extra_sources.tsv for per-protocol sources, and an offline check to
tools/selftest_sources.py for any pure logic (URL/path mapping, parsing).
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from importlib import import_module
from pathlib import Path
from typing import Callable

from .base import ROOT, Folder, say

EXTRA_SOURCES = ROOT / "catalogue" / "extra_sources.tsv"
MODULES = ("scg_github", "protocolsio", "vendor", "elife", "pmc_cloud")


@dataclass
class Handler:
    name: str
    fn: Callable
    doi: re.Pattern | None = None
    url: re.Pattern | None = None
    auto: Callable | None = None
    doc: str = ""


REGISTRY: dict[str, Handler] = {}


def register(name: str, doi: str | None = None, url: str | None = None,
             auto: Callable | None = None, doc: str = ""):
    def deco(fn):
        REGISTRY[name] = Handler(name, fn, re.compile(doi, re.I) if doi else None,
                                 re.compile(url, re.I) if url else None, auto, doc)
        return fn
    return deco


_loaded = False


def load() -> dict[str, Handler]:
    global _loaded
    if not _loaded:
        for m in MODULES:
            import_module(f"{__name__}.{m}")
        _loaded = True
    return REGISTRY


@dataclass
class Job:
    """What a handler gets: the protocol, its folder, and a way back into get_sources."""
    slug: str
    rows: list[dict]
    folder: Folder
    paper: Callable[[str], None] = lambda doi: None
    seen: set = field(default_factory=set)       # DOIs / specs already handled this run


# ------------------------------------------------------------------ extra sources

def parse_extra_sources(text: str) -> list[dict]:
    """catalogue/extra_sources.tsv -> [{slug, source, what, why}]. Raises ValueError on a
    malformed row (so the self-test catches it)."""
    out = []
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        return out
    header = lines[0].split("\t")
    if header[:4] != ["slug", "source", "what", "why"]:
        raise ValueError(f"extra_sources header must be slug/source/what/why, got {header}")
    for i, ln in enumerate(lines[1:], 2):
        f = ln.split("\t")
        if len(f) < 2 or not f[0].strip() or not f[1].strip():
            raise ValueError(f"extra_sources row {i}: needs slug and source: {ln!r}")
        f = (f + ["", ""])[:4]
        out.append(dict(zip(("slug", "source", "what", "why"), (c.strip() for c in f))))
    return out


def extra_sources(slug: str | None = None, path: Path = EXTRA_SOURCES) -> list[dict]:
    if not path.exists():
        return []
    rows = parse_extra_sources(path.read_text(encoding="utf-8"))
    return [r for r in rows if slug is None or r["slug"] == slug]


def split_spec(source: str) -> tuple[str, str]:
    """'scg:data/BD/*' -> ('scg', 'data/BD/*'); a URL -> ('', url)."""
    m = re.match(r"^([a-z][\w-]*):(?!//)(.*)$", source.strip())
    return (m.group(1), m.group(2)) if m else ("", source.strip())


def handler_for_url(url: str) -> Handler | None:
    for h in load().values():
        if h.url and h.url.search(url):
            return h
    return REGISTRY.get("url")


def handlers_for_doi(doi: str) -> list[Handler]:
    return [h for h in load().values() if h.doi and h.doi.search(doi)]


def run_spec(job: Job, source: str, what: str = "", why: str = "") -> None:
    name, arg = split_spec(source)
    key = f"{name}:{arg}"
    if key in job.seen:
        return
    job.seen.add(key)
    desc = "; ".join(x for x in (what, why and f"why: {why}") if x)
    if name in ("doi", "pmc"):
        job.paper(arg)
        return
    if name:
        h = load().get(name)
        if not h:
            say(f"    !! extra source {source!r}: no fetcher named {name!r}")
            return
        h.fn(job, arg, desc)
        return
    h = handler_for_url(arg)
    if h:
        h.fn(job, arg, desc)


def run_auto(job: Job) -> None:
    for h in load().values():
        if h.auto and h.auto(job):
            h.fn(job, "", "")


def run_extra(job: Job) -> None:
    for r in extra_sources(job.slug):
        say(f"  extra  {r['source']}")
        run_spec(job, r["source"], r["what"], r["why"])
