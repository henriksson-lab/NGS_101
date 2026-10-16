"""Check source-exact oligo inputs against locally extracted source documents.

Only graph input materials with an exact sequence state are candidates. A missing
source directory is a SKIP; a present source bundle with no matching sequence is
a failure. This is a transcription guard, not a substitute for reading the source.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))
sys.path.insert(0, str(ROOT / "tools"))

from chemdraw import Construct, Segment  # noqa: E402
from protocol_graph import MaterialKind, ProtocolGraph  # noqa: E402
import scrape_primers  # noqa: E402

BASES = re.compile(r"[ACGTUNRYSWKMBDHV]{12,}", re.I)


@dataclass(frozen=True)
class Candidate:
    material: str
    sequence: str


def _sequences(state) -> tuple[str, ...]:
    """Return whole, ordered, source-exact strands; never construct fragments."""
    if isinstance(state, str):
        return (state,) if BASES.fullmatch(state) else ()
    if isinstance(state, Segment):
        return () if state.placeholder or state.inferred else _sequences(state.top)
    if isinstance(state, Construct):
        parts = state.segments
        if any(p.placeholder or p.inferred for p in parts):
            return ()
        return _sequences("".join(p.top for p in parts))
    if isinstance(state, (tuple, list)):
        if all(isinstance(x, Segment) for x in state):
            if any(x.placeholder or x.inferred for x in state):
                return ()
            return _sequences("".join(x.top for x in state))
        if all(isinstance(x, str) for x in state):
            return tuple(s for x in state for s in _sequences(x))
    return ()


def candidates(graph: ProtocolGraph) -> tuple[Candidate, ...]:
    return tuple(Candidate(m.id, seq) for m in graph.materials.values()
                 if m.id not in graph._producer
                 and m.kind in {MaterialKind.OLIGO, MaterialKind.ADAPTER}
                 and not m.inferred
                 for seq in _sequences(m.state))


def unrepresented_input_ids(graph: ProtocolGraph) -> tuple[str, ...]:
    """Oligo/adapter input IDs without an exact whole-strand state to compare."""
    return tuple(m.id for m in graph.materials.values()
                 if m.id not in graph._producer
                 and m.kind in {MaterialKind.OLIGO, MaterialKind.ADAPTER}
                 and not m.inferred and not _sequences(m.state))


def unrepresented_inputs(graph: ProtocolGraph) -> int:
    return len(unrepresented_input_ids(graph))


def declared_locators(slug: str) -> tuple[str, ...]:
    """Locator fragments for the sources extra_sources.tsv declares for this protocol.

    A declared predecessor, vendor manual or protocol page carries the defining
    oligos for many protocols, and its filename never contains the protocol's own
    DOI, so the DOI filter below would discard exactly the file that was fetched
    to supply them.
    """
    path = ROOT / "catalogue" / "extra_sources.tsv"
    if not path.is_file():
        return ()
    import csv
    result = []
    with path.open(encoding="utf-8", newline="") as handle:
        rows = (line for line in handle if not line.startswith("#"))
        for row in csv.DictReader(rows, delimiter="\t"):
            if row.get("slug") != slug:
                continue
            source = (row.get("source") or "").strip()
            if not source:
                continue
            arg = source.split(":", 1)[1] if ":" in source and "//" not in source else source
            # protocols.io references often end in /v1, while fetched PDFs use
            # the protocol accession in their filename rather than the DOI.
            accession = re.search(r"protocols\.io\.([a-z0-9]+)", arg, re.I)
            for token in (arg, arg.rsplit("/", 1)[-1],
                          accession.group(1) if accession else ""):
                token = token.strip().strip("*").lower()
                if len(token) >= 6:
                    result.append(token)
    return tuple(dict.fromkeys(result))


def source_texts(folder: Path, doi: str = "") -> tuple[tuple[Path, str], ...]:
    """Use the doctext twins and primary plain-text sources; exclude derivatives."""
    if not folder.is_dir():
        return ()
    declared = declared_locators(folder.name)
    # The scg_lib_structs drawn page is inspiration, never evidence. Its content arrives
    # under two prefixes: upstream_* (the page) and scg_docs_* / scg_adapters_* (the same
    # page's markdown and scg's own derived adapter FASTAs). A prefix is not provenance --
    # genuine vendor PDFs also land as upstream_* when the page links them directly -- so
    # the manifest's `what` column is the authority when it matters.
    excluded = ("upstream_", "scg_docs_", "scg_adapters_")
    paths = sorted(p for p in folder.rglob("*.txt")
                   if "_invalid" not in p.parts
                   and not p.name.startswith(excluded)
                   and p.name != "MANIFEST.tsv")
    manifest = folder / "MANIFEST.tsv"
    manifest_urls = {}
    if manifest.exists():
        import csv
        with manifest.open(encoding="utf-8", newline="") as handle:
            manifest_urls = {row["file"]: row["url"] for row in csv.DictReader(handle, delimiter="\t")}
    result = []
    suffix = doi.rsplit("/", 1)[-1].lower() if doi else ""
    for path in paths:
        data = path.read_text(encoding="utf-8", errors="replace")
        original = path.name[:-4]
        locator = (original + " " + manifest_urls.get(original, "")).lower()
        if suffix and suffix not in locator and doi.lower() not in data[:5000].lower() \
                and not any(token in locator for token in declared):
            continue
        result.append((path, data))
    return tuple(result)


def excluded_text_count(folder: Path, doi: str = "") -> int:
    """How many extracted text files the locator filter dropped, for diagnostics.

    A predecessor or vendor document that is physically present but undeclared is
    discarded by the filter, which makes "the oligo is in a file we refused to read"
    indistinguishable from "no source prints this oligo". Reporting the count turns the
    second reading into a checkable first one.
    """
    if not folder.is_dir():
        return 0
    present = [p for p in folder.rglob("*.txt")
               if "_invalid" not in p.parts
               and not p.name.startswith(("upstream_", "scg_docs_", "scg_adapters_"))
               and p.name != "MANIFEST.tsv"]
    return max(0, len(present) - len(source_texts(folder, doi)))


def _exact_hit(text: str, seq: str, scan=None) -> bool:
    """Require every fixed query base to agree; N in the model is variable."""
    for _, strand, _, written in scrape_primers.find(text, seq, sc=scan, detail=True):
        try:
            observed = scrape_primers._letters(written)
        except ValueError:
            continue
        if strand == "rc":
            observed = scrape_primers.revcomp(observed)
        elif strand == "rev":
            observed = observed[::-1]
        elif strand == "comp":
            observed = scrape_primers.complement(observed)
        if len(observed) == len(seq) and all(a == "N" or a == b
                                             for a, b in zip(seq.upper().replace("U", "T"), observed)):
            return True
    return False


def check_graph(graph: ProtocolGraph, folder: Path,
                exceptions: dict[str, str] | None = None,
                doi: str = "") -> tuple[list[str], list[str], int]:
    """Return (failures, skips, checked count). Exceptions key graph material IDs."""
    exceptions = exceptions or {}
    found = candidates(graph)
    unknown = set(exceptions) - {c.material for c in found}
    failures = [f"{graph.name}: stale exception for {x}" for x in sorted(unknown)]
    failures += [f"{graph.name}: exception for {x} needs a reason"
                 for x, reason in exceptions.items() if not reason.strip()]
    if not found:
        return failures, [], 0
    docs = source_texts(folder, doi)
    if not docs:
        return failures, [f"{graph.name}: no extracted source text in {folder}"], 0
    scans = [(path, text, scrape_primers.scan_copy(text)) for path, text in docs]
    checked = 0
    for candidate in found:
        if candidate.material in exceptions:
            continue
        checked += 1
        if not any(_exact_hit(text, candidate.sequence, scan)
                   for _, text, scan in scans):
            note = ""
            dropped = excluded_text_count(folder, doi)
            if dropped:
                note = (f" ({dropped} text file(s) in the bundle excluded by the locator"
                        " filter; declare the source in catalogue/extra_sources.tsv if it"
                        " carries this oligo)")
            failures.append(f"{graph.name}: {candidate.material} {candidate.sequence} "
                            f"absent from local source text{note}")
    return failures, [], checked
