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

from chemdraw import Construct, Scene, Segment  # noqa: E402
from protocol_graph import MaterialKind, ProtocolGraph  # noqa: E402
import scrape_primers  # noqa: E402

BASES = re.compile(r"[ACGTUNRYSWKMBDHV]{12,}", re.I)


@dataclass(frozen=True)
class Candidate:
    material: str
    sequence: str


MIN_FIXED = 12   # concrete bases a strand with variable windows must still pin down


def _biorxiv_doi(doi: str) -> bool:
    # Cold Spring Harbor also publishes Genome Research under 10.1101/gr.*.
    return re.match(r"10\.1101/\d", doi.lower()) is not None


_CLASSES = set("NRYSWKMBDHV")


def _window(seg: Segment) -> int | None:
    """Length of a variable placeholder window, or None if it is not one.

    Only an identifier window (a segment with a declared barcode/index/UMI feature) or a
    degenerate IUPAC-class placeholder such as "VN" or "B" is a variable window. Any
    other placeholder (an insert, a template or cDNA stand-in) is not reagent sequence.
    """
    top = seg.top.upper()
    degenerate = bool(top) and set(top) <= _CLASSES and bool(set(top) - {"N"})
    if seg.feature is None and not degenerate:
        return None
    if seg.length_bp:
        return seg.length_bp
    return None if seg.is_role_token() or not seg.top else len(seg.top)


# Modifications that carry no nucleotide: a placeholder written only in these tokens
# (a C3 spacer between two docking arms) adds no base to the strand.
_NON_NUCLEOTIDE = {"spacer", "biotin", "amine", "phosphate", "thiol", "acrydite"}


def _notation_bases(seg: Segment) -> str | None:
    """Bases a placeholder written purely in source modification notation stands for.

    "/rG//rG//iXNA_G/" is three modified G; "[C3]" or "/iSpC3/" is a non-nucleotide
    spacer and adds none. Tokens are read with the scraper's own vocabulary, so the
    model and the printed oligo are compared through one parser. None when the text is
    not wholly such tokens, a token is unknown, or a declared length disagrees.
    """
    top = seg.top.strip()
    if not re.fullmatch(r"(?:/{1,2}[\w-]+)+/?|(?:\[[\w-]+\])+", top):
        return None
    bases = []
    for token in re.findall(r"[^/\[\]]+", top):
        if re.fullmatch(r"[rm+][ACGTU]", token):
            bases.append(token[-1])
        elif scrape_primers.internal_base(token):
            bases.append(scrape_primers.internal_base(token))
        elif scrape_primers.canonical(token) not in _NON_NUCLEOTIDE:
            return None
    seq = "".join(bases)
    if seg.length_bp is not None and seg.length_bp != len(seq):
        return None
    return seq


def _strand(segs) -> tuple[str, ...]:
    """One whole strand; placeholder windows of fixed length become variable N positions.

    A placeholder's letters are a stand-in (AAAAAAAA for "barcode A", "B" or "VN" as an
    IUPAC class), never bases to compare, so every placeholder position is N, which
    _exact_hit treats as variable. An inferred segment or a window of unknown length
    leaves the strand without an exact state.
    """
    parts = []
    windows = False
    for seg in segs:
        if seg.inferred:
            return ()
        if seg.placeholder:
            n = _window(seg)
            if n is None:
                written = _notation_bases(seg)
                if written is None:
                    return ()
                parts.append(written)
                continue
            parts.append("N" * n)
            windows = True
        else:
            parts.append(seg.top)
    seq = "".join(parts)
    if windows and sum(seq.upper().count(b) for b in "ACGTU") < MIN_FIXED:
        return ()
    return _sequences(seq)


def _sequences(state) -> tuple[str, ...]:
    """Return whole, ordered, source-exact strands; never construct fragments."""
    if isinstance(state, str):
        return (state,) if BASES.fullmatch(state) else ()
    if isinstance(state, Segment):
        return _strand((state,))
    if isinstance(state, Construct):
        return _strand(state.segments)
    if isinstance(state, Scene):
        # Scene keeps the original ordering strands separately from their drawn
        # orientation and pairing offsets. Never join different reagent strands
        # or reverse a strand merely because it is drawn antiparallel.
        return tuple(dict.fromkeys(seq for strand in state.strands.values()
                                   for seq in _strand(strand.segs)))
    if isinstance(state, (tuple, list)):
        if all(isinstance(x, Segment) for x in state):
            return _strand(state)
        if all(isinstance(x, str) for x in state):
            return tuple(s for x in state for s in _sequences(x))
        if all(isinstance(x, (Scene, Construct)) for x in state):
            # separate reagents delivered together (a loaded transposome's two adapters)
            return tuple(dict.fromkeys(s for x in state for s in _sequences(x)))
    return ()


def candidates(graph: ProtocolGraph) -> tuple[Candidate, ...]:
    return tuple(Candidate(m.id, seq) for m in graph.materials.values()
                 if m.id not in graph._producer
                 and m.kind in {MaterialKind.OLIGO, MaterialKind.ADAPTER}
                 and not m.inferred
                 for seq in _sequences(m.state))


def unrepresented_input_ids(graph: ProtocolGraph) -> tuple[str, ...]:
    """Oligo/adapter inputs with at least one unrepresented whole strand."""
    def incomplete(state) -> bool:
        if isinstance(state, Scene):
            return not state.strands or any(not _strand(strand.segs)
                                            for strand in state.strands.values())
        if (isinstance(state, (tuple, list)) and state
                and all(isinstance(x, (Scene, Construct)) for x in state)):
            return any(incomplete(x) for x in state)
        return not _sequences(state)

    return tuple(m.id for m in graph.materials.values()
                 if m.id not in graph._producer
                 and m.kind in {MaterialKind.OLIGO, MaterialKind.ADAPTER}
                 and not m.inferred and incomplete(m.state))


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
            # "handler:arg"; the arg may itself be a URL ("protocolsio:https://...")
            handler = re.match(r"[a-z]+:(?!//)", source)
            arg = source[handler.end():] if handler else source
            # protocols.io references often end in /v1, while fetched PDFs use
            # the protocol accession in their filename rather than the DOI.
            accession = re.search(r"protocols\.io\.([a-z0-9]+)", arg, re.I)
            for token in (arg, arg.rstrip("/").rsplit("/", 1)[-1],
                          accession.group(1) if accession else ""):
                token = token.strip().strip("*").lower()
                if len(token) >= 6:
                    result.append(token)
    return tuple(dict.fromkeys(result))


def linked_version_dois(slug: str, manifest: Path | None = None) -> tuple[str, ...]:
    """Other DOIs of this protocol's defining work: catalogue `primary` papers of its directory.

    catalogue/scg_lib_structs.tsv lists a protocol's published article and its preprint
    as primary papers of one directory, and get_sources.py fetches such a paper only
    after linking it to the defining DOI; a paper it could not link gets an
    "(note) unlinked catalogue paper DOI" manifest row and is left out here. Without
    this a directory named for a preprint would never read its published version.
    """
    path = ROOT / "catalogue" / "scg_lib_structs.tsv"
    if not path.is_file():
        return ()
    unlinked = ""
    if manifest is not None and manifest.is_file():
        unlinked = manifest.read_text(encoding="utf-8", errors="replace").lower()
    import csv
    found = []
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            doi = (row.get("doi") or "").strip().lower()
            if (row.get("our_dir") == slug and row.get("role") == "primary" and doi
                    and f"(note) unlinked catalogue paper {doi}" not in unlinked):
                found.append(doi)
    # A linked preprint is historical once this defining work has a published DOI.
    # Its oligos may only be handled through an explicitly justified local exception,
    # not silently accepted as a transcription match for the published protocol.
    published = [d for d in found if not _biorxiv_doi(d)]
    return tuple(dict.fromkeys(published or found))


def manifest_locators(manifest: Path) -> dict[str, tuple[str, ...]]:
    """File name -> its manifest URL(s), plus the name and URL of byte-identical copies.

    get_sources.py records a second download of the same bytes as a "(duplicate) NAME"
    row whose `what` ends "identical to TWIN", and keeps only TWIN's text. A published
    supplement identical to its preprint copy, or a publisher file identical to a
    mirrored copy, is thereby identified by the duplicate row's name and URL; both
    files of such a pair carry the other's locator.
    """
    if not manifest.is_file():
        return {}
    import csv
    found: dict[str, list[str]] = {}
    pairs = []
    with manifest.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            name = (row.get("file") or "").strip()
            duplicate = name.startswith("(duplicate)")
            name = name[len("(duplicate)"):].strip() if duplicate else name
            url = (row.get("url") or "").strip()
            found.setdefault(name, []).append(url)
            twin = re.search(r"identical to (\S+)\s*$", row.get("what") or "")
            if duplicate and twin:
                pairs.append((name, twin.group(1)))
    for name, twin in pairs:
        own, other = list(found.get(name, ())), list(found.get(twin, ()))
        found.setdefault(twin, []).extend([name, *own])
        found[name].extend([twin, *other])
    return {name: tuple(dict.fromkeys(x for x in urls if x)) for name, urls in found.items()}


# The scg_lib_structs drawn page is inspiration, never evidence. Its content arrives
# under two prefixes: upstream_* (the page) and scg_docs_* / scg_adapters_* (the same
# page's markdown and scg's own derived adapter FASTAs). A prefix is not provenance --
# genuine vendor PDFs also land as upstream_* when the page links them directly -- so
# the manifest's `what` column is the authority when it matters.
EXCLUDED = ("upstream_", "scg_docs_", "scg_adapters_")
# A FASTA original (a vendor primer file) is plain text with no doctext twin. FASTA files
# in the scg_lib_structs repository are scg's own derived tables, so scg_* FASTAs stay out.
FASTA = {".fa", ".fasta", ".fna"}


def _present(folder: Path) -> list[Path]:
    """Candidate source texts: doctext twins and primary FASTA files, derivatives excluded."""
    texts = [p for p in folder.rglob("*.txt") if p.name != "MANIFEST.tsv"]
    have = {p.name for p in texts}
    fasta = [p for p in folder.rglob("*") if p.suffix.lower() in FASTA and p.is_file()
             and p.name + ".txt" not in have and not p.name.startswith("scg_")]
    return sorted(p for p in texts + fasta
                  if "_invalid" not in p.parts and not p.name.startswith(EXCLUDED))


def source_texts(folder: Path, doi: str = "") -> tuple[tuple[Path, str], ...]:
    """Use the doctext twins and primary plain-text sources; exclude derivatives."""
    if not folder.is_dir():
        return ()
    declared = declared_locators(folder.name)
    paths = _present(folder)
    aliases = manifest_locators(folder / "MANIFEST.tsv")
    result = []
    dois = [doi.lower()] if doi else []
    if doi:
        dois += [d for d in linked_version_dois(folder.name, folder / "MANIFEST.tsv")
                 if d not in dois]
        published = [d for d in dois if not _biorxiv_doi(d)]
        dois = published or dois
    suffixes = [d.rsplit("/", 1)[-1] for d in dois]
    suffixes = suffixes[:1] + [x for x in suffixes[1:] if len(x) >= 6]
    for path in paths:
        data = path.read_text(encoding="utf-8", errors="replace")
        original = path.name[:-4] if path.suffix == ".txt" else path.name
        # a primary plain-text original (a GEO SOFT text) is its own manifest entry
        names = dict.fromkeys((original, path.name))
        locator = " ".join([*names, *(a for n in names for a in aliases.get(n, ()))]).lower()
        head = data[:5000].lower()
        if dois and not any(x in locator for x in suffixes) \
                and not any(d in head for d in dois) \
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
    return max(0, len(_present(folder)) - len(source_texts(folder, doi)))


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
        if _agrees(seq, observed):
            return True
    return False


def _agrees(seq: str, observed: str) -> bool:
    """Every fixed model base equals the printed letter; a model N is variable.

    A printed placeholder of unstated length ("[i7]", read as "…") may stand only for one
    run of model N positions: the fixed bases on both sides must still agree in full.
    """
    seq = seq.upper().replace("U", "T")
    if scrape_primers.GAP not in observed:
        return len(observed) == len(seq) and all(a == "N" or a == b
                                                 for a, b in zip(seq, observed))
    pattern = "".join("N+" if c == scrape_primers.GAP else f"[{re.escape(c)}N]"
                      for c in observed)
    return re.fullmatch(pattern, seq) is not None


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
