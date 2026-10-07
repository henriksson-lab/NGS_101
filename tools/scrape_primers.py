#!/usr/bin/env python3
"""Pull every oligo-looking sequence, with the sentence around it, out of source documents.

A triage aid for reading a paper, a supplement or a scraped method page: instead of
paging through 40 pages for the oligo table, list each place a DNA sequence appears,
what it is called nearby, and which sequences we already model it contains. Heuristic
by design -- a hit is a pointer into the document, never a fact. Verify against the
source before anything goes into a protocol module.

    python3 tools/scrape_primers.py pdf/hagemann.txt
    python3 tools/scrape_primers.py DIR                     # every readable file under it
    python3 tools/scrape_primers.py doc.pdf --tsv > hits.tsv
    python3 tools/scrape_primers.py doc.pdf --known-only    # only hits containing a lib sequence
    python3 tools/scrape_primers.py DIR --all               # no family collapsing, no --max-hits cut
    python3 tools/scrape_primers.py DIR --find AGATGTGTATAAGAGACAG   # where is this oligo?
    python3 tools/scrape_primers.py DIR --find "5'-/5Phos/CAGAGCNNNNNNNN[10bp barcode]T30-3'"

Reads .pdf (via poppler's `pdftotext`), .docx/.xlsx (standard library), .html/.htm, and
plain text (.txt .md .tsv .csv .json .fa/.fasta .fq/.fastq). In a directory a document
with a `FILE.txt` twin (tools/doctext.py) is read as the twin, so line numbers point into
it; a file with the same content as one already read is skipped. A .pdf/.xlsx that is
really an HTML page (a captcha or login page saved by a fetcher) is skipped with that
reason.

What counts as one sequence (the scan reads a normalised copy of the text and maps every
hit back to the original, which is reported "as written"):
  * IUPAC bases, upper or lower case: CAAGCAGAAGACGGCATACGAGATccgaatccgaGTCTCGTGGGCTCGG
    is one oligo. A run with no upper-case base counts only as a whole table cell
    (prose like "acgtacgt" in a sentence does not).
  * pieces rejoined across a PDF line wrap, a stray space or hyphen, codon-spaced
    triplets ("TCA GAC GTG TGC ..."), and a PDF line number between the halves. Not
    rejoined: lines indented under one another (a drawn duplex), and runs of 3 or more
    lines that are each one sequence of the same length (a barcode list).
  * per-base marks are looked through: RNA rN, 2'-OMe mN, LNA +N, phosphorothioate * or
    ∗ (also "-s-"), IDT internal /iXxx/ (/ideoxyU/ reads as U), (dU), iso-bases iCiGiC.
  * shorthand and placeholders inside an oligo keep it whole: T30, T30VN, dT30, (T)30,
    (dT)30, T(30), N8, N16-N12, (N)8, XXXXX and JJJJJ (one N each), [NNNNNN],
    [8-bp Round1 barcode], [10bp barcode], [UMI10], [BC6], [T30], {UMI}, <cbc> are
    expanded to that many bases in `sequence`; a placeholder without a length ([i7],
    [barcode], XXX...XXX, NNNN...NNNN, (dT), (A)n, cDNA_Insert) becomes "…" in
    `sequence` and the hit's length is shown as "N+ nt".

Each hit lists: the bases; the oligo as written with its modifications; the parsed
modifications, named the same way however they are written ("5' biotin (/5Biosg/)",
"5' biotin (Bio)", "5' phosphate ([5Phos])", "3' rGrGrG", "RNA (written with U)", and
"written 3'->5'" for a strand drawn 3'-...-5'); a nearby name (in a table, the row's
identifier cell; else the token before it or a short line above it); the lib/ sequences
it contains on either strand (whole, or a >= 15-nt piece "(nt a-b of L)"); other lines
with the same oligo; and the sentence around it (for a table row, the row). Tuned for
recall: expect false hits.

Large tables: 4+ consecutive rows that are one oligo with a varying barcode or index
(same length and differing at a minority of positions, or the same name stem -- P7-1,
P7-2, ... -- and a shared 12-nt stretch) are collapsed into one hit, "family: 96
variants, lines a-b, variable at nt i-j" -- `--all` lists every row. Files are ordered by
their most oligo-like hit and each hit has a `score`; `--max-hits` (default 200 per file,
counted after collapsing) keeps the highest scorers and says how many it cut. A 30 MB
table text takes a few seconds, a 100 MB guide library about a minute.

--find SEQ (repeatable) lists where SEQ occurs. The query may be pasted as written:
5'/3' ends, /mods/, rN/mN/+N/* marks, spaces, hyphens, shorthand and placeholders are
read as above. Letters are base classes on both sides, case-insensitive: N matches any
base, U = T, R = A/G, ...; a query N-run also matches a "…" placeholder in the text
(--strict: letters must be equal, so N matches only N). Line wraps, spaces, hyphens,
tabs and modification marks in the text are skipped. Each location is tagged with the
strand it was found on (+, rc; rev and comp = a strand written 3'->5') and "exact" or
"iupac" (matched only through a degenerate letter, e.g. a query N on a concrete base);
a place where under 60% of the query's A/C/G/T meet the same letter (a drawn NNN...N) is
not listed. Per file, exact locations come first; --max-hits caps the lines shown.
"""
from __future__ import annotations

import argparse
import hashlib
import html
import importlib
import re
import shutil
import subprocess
import sys
from bisect import bisect_left, bisect_right
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))

TEXT_EXT = {".txt", ".md", ".tsv", ".csv", ".json", ".fa", ".fasta", ".fq", ".fastq"}
OFFICE_EXT = {".xlsx", ".xlsm", ".docx"}
HTML_EXT = {".html", ".htm"}

BASE = "ACGTUNRYSWKMBDHV"
LOW = "acgtun"
GAP = "…"           # an uncounted placeholder, in the scan copy and in Hit.seq
HARD_NL = "\v"      # a line break the scan must not join across (scan copy only)

# Thresholds lean towards recall: the output is read by a person or an LLM who can discard
# a false hit in a second, whereas a missed oligo is never seen at all.
MIN_UPPER = 6
MIN_LEN = 12
MIN_CONCRETE = 8         # A/C/G/T/U bases; a run of N alone is not an oligo
MIN_ACGT = 0.6           # of the letters outside N placeholders: motifs (TATGCWAATBAV) are not

KEYWORDS = re.compile(r"primer|oligo|adapt(?:e|o)r|probe|barcode|index|linker|TSO|"
                      r"template.switch|RT\b|PCR|read ?[12]|i[57]\b|P[57]\b|UMI", re.I)

_IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "U": "T", "R": "AG", "Y": "CT", "S": "CG",
          "W": "AT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG",
          "N": "ACGT"}
_COMP = str.maketrans("ACGTURYSWKMBDHVNacgturyswkmbdhvn", "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn")


def complement(seq: str) -> str:
    return seq.translate(_COMP)


def revcomp(seq: str) -> str:
    """Reverse complement; IUPAC letters complement to their IUPAC partner."""
    return complement(seq)[::-1]


# --------------------------------------------------------------------- the hits

@dataclass
class Hit:
    path: Path
    line: int
    seq: str                              # bases as the scan read them (case kept, "…" = gap)
    sentence: str
    label: str = ""
    known: list[str] = field(default_factory=list)
    written: str = ""                     # as in the source, modifications included
    mods: list[str] = field(default_factory=list)
    also: list[int] = field(default_factory=list)    # further lines with the same oligo
    table_cell: bool = False
    family: list["Hit"] = field(default_factory=list)  # collapsed rows, this one first
    windows: list[tuple[int, int]] = field(default_factory=list)  # 1-based nt that vary

    @property
    def keyword(self) -> bool:
        return bool(KEYWORDS.search(self.sentence))

    @property
    def length(self) -> int:
        return len(self.seq.replace(GAP, ""))

    @property
    def length_text(self) -> str:
        return f"{self.length}{'+' if GAP in self.seq else ''} nt"

    @property
    def score(self) -> int:
        """How oligo-like the hit is, for ordering and the --max-hits cut.
        Up: modifications, a sequence we already model, primer vocabulary nearby, a name.
        Down: degenerate IUPAC letters, which motif tables (TATGCWAATBAV) are full of and
        real oligos almost never are -- an oligo-dT's trailing VN aside."""
        s = 3 * bool(self.mods) + 3 * bool(self.known) + 2 * self.keyword + bool(self.label)
        s += 18 <= self.length <= 200
        bare = re.sub(r"([RYSWKMBDHV])\1{3,}", "", self.seq.rstrip("VNvn"))  # YYYYYYYY = an index
        degenerate = sum(c in "RYSWKMBDHVrysw" for c in bare)
        return s - 2 * min(degenerate, 3) + bool(self.family)


# ---------------------------------------------------------------------- reading

def html_text(raw: str) -> str:
    raw = re.sub(r"(?is)<(script|style)\b.*?</\1>", " ", raw)
    # block boundaries become line breaks; inline tags (the coloured <span>s a sequence is
    # split into on scg_lib_structs pages) vanish, so the bases rejoin
    raw = re.sub(r"(?i)<br\s*/?>|</(?:p|div|tr|li|h[1-6]|pre|td|th|seq|align|table)>", "\n", raw)
    raw = re.sub(r"<[^>]+>", "", raw)
    return html.unescape(raw)


def _xml_text(xml: str) -> str:
    return html.unescape(re.sub(r"<[^>]+>", "", xml))


def xlsx_text(path: Path) -> str:
    """Every sheet as tab-separated rows, under a '## sheet' heading. Standard library only."""
    import zipfile
    with zipfile.ZipFile(path) as z:
        names = z.namelist()
        shared = []
        if "xl/sharedStrings.xml" in names:
            sst = z.read("xl/sharedStrings.xml").decode("utf-8", "replace")
            shared = [_xml_text(si) for si in re.findall(r"<si>(.*?)</si>", sst, re.S)]
        book = z.read("xl/workbook.xml").decode("utf-8", "replace") if "xl/workbook.xml" in names else ""
        titles = re.findall(r'<sheet [^>]*name="([^"]*)"', book)
        sheets = sorted((n for n in names if re.match(r"xl/worksheets/sheet\d+\.xml$", n)),
                        key=lambda n: int(re.search(r"(\d+)", n.rsplit("/", 1)[1]).group(1)))
        out = []
        for i, n in enumerate(sheets):
            out.append(f"## sheet {html.unescape(titles[i]) if i < len(titles) else i + 1}")
            xml = z.read(n).decode("utf-8", "replace")
            for row in re.findall(r"<row\b[^>]*>(.*?)</row>", xml, re.S):
                cells = []
                for attrs, body in re.findall(r"<c\b([^>]*?)(?:/>|>(.*?)</c>)", row, re.S):
                    v = re.search(r"<v>(.*?)</v>", body or "", re.S)
                    if 't="s"' in attrs and v:
                        idx = int(v.group(1))
                        cells.append(shared[idx] if idx < len(shared) else "")
                    elif 't="inlineStr"' in attrs:
                        cells.append(_xml_text(body))
                    else:
                        cells.append(html.unescape(v.group(1)) if v else "")
                if any(c.strip() for c in cells):
                    out.append("\t".join(c.replace("\n", " ") for c in cells))
        return "\n".join(out) + "\n"


def docx_text(path: Path) -> str:
    """Paragraphs one per line; table cells tab-separated."""
    import zipfile
    with zipfile.ZipFile(path) as z:
        xml = z.read("word/document.xml").decode("utf-8", "replace")
    # inside a table cell, paragraphs join with a space so one row stays one line
    xml = re.sub(r"<w:tc\b.*?</w:tc>",
                 lambda m: re.sub(r"</w:p>|<w:br/>", " ", m.group(0)).replace("</w:tc>", "\t"),
                 xml, flags=re.S)
    xml = re.sub(r"</w:p>|</w:tr>|<w:br/>", "\n", xml)
    xml = re.sub(r"<w:tab/>", "\t", xml)
    return _xml_text(xml)


def _looks_html(head: bytes) -> bool:
    h = head.lstrip().lower()
    return h.startswith((b"<!doctype html", b"<html")) or (h.startswith(b"<") and b"<html" in h)


def read(path: Path) -> str | None:
    """The document as plain text, or None (with the reason on stderr) if it cannot be read."""
    ext = path.suffix.lower()
    if ext in OFFICE_EXT | {".pdf"}:
        with open(path, "rb") as fh:
            head = fh.read(4096)
        if _looks_html(head):
            print(f"skip {path}: an HTML page saved as {ext}, not a {ext[1:].upper()} "
                  f"(a captcha, login or error page?) -- fetch the file by hand", file=sys.stderr)
            return None
    if ext in OFFICE_EXT:
        try:
            return xlsx_text(path) if ext != ".docx" else docx_text(path)
        except Exception as e:                          # noqa: BLE001
            print(f"skip {path}: cannot read ({e})", file=sys.stderr)
            return None
    if ext == ".pdf":
        return _pdf(path)
    raw = path.read_bytes()
    if raw[:5] == b"%PDF-":                 # a PDF saved under another name
        return _pdf(path)
    if b"\0" in raw[:4096]:
        print(f"skip {path}: binary, not a format this reads", file=sys.stderr)
        return None
    text = raw.decode("utf-8", errors="replace")
    if ext in HTML_EXT or text.lstrip()[:15].lower().startswith(("<!doctype", "<html")):
        return html_text(text)
    return text


def _pdf(path: Path) -> str | None:
    if not shutil.which("pdftotext"):
        print(f"skip {path}: pdftotext (poppler-utils) not installed", file=sys.stderr)
        return None
    done = subprocess.run(["pdftotext", str(path), "-"], capture_output=True, text=True)
    if done.returncode:
        print(f"skip {path}: pdftotext failed", file=sys.stderr)
        return None
    return done.stdout


def files(args: list[str]) -> list[Path]:
    """Files to read. In a directory, a document with a converted twin (paper.pdf beside
    paper.pdf.txt, from tools/doctext.py) is read as the text, so hits cite its lines."""
    out = []
    for a in args:
        p = Path(a)
        if p.is_dir():
            found = sorted(f for f in p.rglob("*") if f.is_file() and
                           f.suffix.lower() in TEXT_EXT | HTML_EXT | OFFICE_EXT | {".pdf"})
            have = {f.name for f in found}
            out += [f for f in found if f.name + ".txt" not in have]
        elif p.exists():
            out.append(p)
        else:
            print(f"skip {a}: no such file", file=sys.stderr)
    return out


# --------------------------------------------------------------- known sequences

def known_sequences() -> dict[str, str]:
    """NAME -> sequence for every upper-case DNA constant (>= 12 nt) in lib/."""
    out: dict[str, str] = {}
    for f in sorted((ROOT / "lib").glob("*.py")):
        try:
            mod = importlib.import_module(f.stem)
        except Exception:                                  # noqa: BLE001
            continue
        for k, v in vars(mod).items():
            if k.isupper() and isinstance(v, str) and re.fullmatch(r"[ACGT]{12,}", v):
                out.setdefault(v, f"{f.stem}.{k}")
    return {name: seq for seq, name in out.items()}


_KMER = 15
_NOT_ACGT = re.compile(f"[^ACGT{GAP}]")
_index_cache: dict[int, dict[str, list[tuple[str, int]]]] = {}


def _kmer_index(known: dict[str, str]) -> dict[str, list[tuple[str, int]]]:
    idx = _index_cache.get(id(known))
    if idx is None:
        idx = {}
        for name, k in known.items():
            for i in range(len(k) - _KMER + 1):
                idx.setdefault(k[i:i + _KMER], []).append((name, i))
        _index_cache.clear()
        _index_cache[id(known)] = idx
    return idx


def contained(seq: str, known: dict[str, str], min_overlap: int = _KMER) -> list[str]:
    """Known sequences in `seq`, either strand: whole ("illumina.P5"), or a shared stretch
    of >= min_overlap nt ("illumina.TRUSEQ_READ2 (nt 15-34 of 34)") -- a primer is often a
    fragment of a modelled sequence, or carries one only in part."""
    s = _NOT_ACGT.sub("N", seq.upper().replace("U", "T"))
    if len(s) < 12:
        return []
    idx = _kmer_index(known)
    found = []                              # (strand, name, a, b, ka, kb) on that strand
    short = [(n, k) for n, k in known.items() if len(k) < min_overlap]
    for strand, x in (("", s), (" (rc)", revcomp(s).replace(GAP, "N"))):
        full = []
        for name, k in short:
            at = x.find(k)
            if at >= 0:
                full.append((at, at + len(k)))
                found.append((strand, name, at, at + len(k), 0, len(k)))
        best: dict[tuple[str, int], tuple[int, int, int, int]] = {}
        for i in range(len(x) - min_overlap + 1):
            for name, off in idx.get(x[i:i + min_overlap], ()):
                d = i - off
                k = known[name]
                if (name, d) in best:
                    continue
                a, ka = i, off              # extend the diagonal both ways
                while a > 0 and ka > 0 and x[a - 1] == k[ka - 1]:
                    a, ka = a - 1, ka - 1
                b, kb = i + min_overlap, off + min_overlap
                while b < len(x) and kb < len(k) and x[b] == k[kb]:
                    b, kb = b + 1, kb + 1
                best[(name, d)] = (a, b, ka, kb)
        for (name, _), (a, b, ka, kb) in best.items():
            if kb - ka == len(known[name]):
                full.append((a, b))
        for (name, _), (a, b, ka, kb) in best.items():
            if kb - ka < len(known[name]) and any(fa <= a and b <= fb for fa, fb in full):
                continue                    # inside a whole match of another sequence
            found.append((strand, name, a, b, ka, kb))
    # lib often holds both strands (P5 and P5_RC): one strand's match is enough. Of
    # overlapping pieces keep the longest; a short whole constant (STEM, 12 nt) inside a
    # longer match adds nothing.
    L = len(s)
    rows = []
    for st, n, a, b, ka, kb in found:
        whole = kb - ka == len(known[n])
        fa, fb = (a, b) if not st else (L - b, L - a)
        cat = 1 if not whole else 0 if len(known[n]) >= _KMER else 2
        rows.append((cat, -(b - a), st, len(known[n]), n, fa, fb, ka, kb, whole))
    rows.sort()
    taken: list[tuple[int, int]] = []
    names: set[str] = set()
    out = []
    for _, _, st, _, n, fa, fb, ka, kb, whole in rows:
        k = known[n]
        if n in names or any(known.get(m) == revcomp(k) for m in names):
            continue
        inside = any(min(fb, b2) - max(fa, a2) >= (0.8 * (fb - fa) if not whole else fb - fa)
                     for a2, b2 in taken)
        if inside and (not whole or len(k) < _KMER):
            continue
        names.add(n)
        taken.append((fa, fb))
        out.append((not whole, bool(st), fa,
                    f"{n}{st}" if whole else f"{n}{st} (nt {ka + 1}-{kb} of {len(k)})"))
    return [t[-1] for t in sorted(out)]


# ---------------------------------------------------------------- modifications
#
# Bases are found on a copy of the text with the per-base marks taken out and the
# shorthand expanded, so that "rGrUrUrCrArG", "A*C*G" or "GTCGACT30VN" reads as one run;
# every hit is then mapped back to the original text, where the oligo is reported exactly
# as written, modifications included.

# Per-base marks: RNA "rG", 2'-OMe "mA", LNA "+A", phosphorothioate "A*C". Four marked
# bases in a row, so prose ("rA", "mM") is left alone.
_MARKED = rf"[rm+][{BASE}][*∗]?(?:[ \t]*\n[ \t]*)?"
PER_BASE = re.compile(rf"{_MARKED}(?:{_MARKED}){{3,}}")
IDT_INTERNAL = re.compile(rf"(?<=[{BASE}{LOW}\]])/(i?[A-Za-z][\w-]*)/(?=[{BASE}{LOW}\[])")

# A 5' modification written directly in front of the bases, and a 3' one directly after.
_MODNAME = (r"/[5-9i]?[\w-]+/|rApp|App|NH2|Amino(?:\s?C\d+|-?linker)?|Am(?:MC6|C6)|"
            r"[Bb]iotin(?:-?TEG|ylated)?|Bio(?:sg|TEG)?|Btn|[Pp]hosphat(?:e|ed)|[Pp]hosphorylated|Phos|Pho|P|FAM|HEX|ROX|TAMRA|Cy[35](?:\.5)?|"
            r"DBCO|Acrydite|Acryd|Azide|PC[ -]?[Ss]pacer|Spacer\s?C\d+|SpC3|C3|C6|ddC|ddN|"
            r"dd[ACGT]|InvdT|idT|BHQ-?[12]|SH|Thiol|OH|(?:i[CG]){2,}")
_MOD = rf"[\[(](?:[53]['′’]?[ -]?)?(?:{_MODNAME})[\])]|{_MODNAME}|(?:[rm+][ACGUTN][*∗]?)+|[*∗]"
# a modification starts a token: "TSO-P AAGC..." is a name, not a 5' phosphate
MOD5 = re.compile(rf"(?:5['′’][ \t]*-?[ \t]*)?(?:(?<![A-Za-z0-9])(?<![A-Za-z0-9]-)(?:{_MOD})"
                  rf"[ \t]*-?[ \t]*)+$")
MOD3 = re.compile(rf"^(?:[ \t]*-?[ \t]*(?:{_MOD})(?![A-Za-z0-9]))+")
_ANY_MARK = re.compile(r"[rm+*∗(/]|-s-")
_R_MARK = re.compile(rf"(?<![A-Za-z])r[{BASE}]|(?<=[{BASE}])r[{BASE}]")
_M_MARK = re.compile(rf"(?<![A-Za-z])m[{BASE}]|(?<=[{BASE}])m[{BASE}]")
_L_MARK = re.compile(rf"\+[{BASE}]")
_TAG = {"r": "RNA bases", "m": "2'-OMe", "+": "LNA", "*": "phosphorothioate"}

# One name per modification, however the paper writes it.
_CANON = [(r"bio|biotin|btn|biosg|bioteg|biotinteg|biotin-?teg|biotinylated", "biotin"),
          (r"phos|pho|p|phosphate|phosphated|phosphorylated", "phosphate"), (r"r?app", "adenylated"),
          (r"nh2|amino.*|ammc6|amc6|aminolinker", "amine"),
          (r"acryd|acrydite", "acrydite"), (r"dd[acgtn]", "dideoxy"),
          (r"sp\d+|spc3|c3|c6|pc[ -]?spacer|spacer\s?c\d+", "spacer"),
          (r"invdt|idt", "inverted dT"), (r"sh|thiol|thiomc6-d|dithiol", "thiol"),
          (r"(?:i[cg]){2,}", "iso-dC/iso-dG")]


_CANON_RX = [(re.compile(rx), name) for rx, name in _CANON]
_BRACKETS = re.compile(r"^[\[(/]+|[\])/]+$")
_END_PREFIX = re.compile(r"^[35]['′’]?|^i(?=[A-Z])")


def canonical(mod: str) -> str:
    """'/5Biosg/' -> 'biotin', '[5Phos]' -> 'phosphate', 'rApp' -> 'adenylated'; else ''."""
    key = _BRACKETS.sub("", mod.strip())
    if re.fullmatch(r"(?:i[CG]){2,}", key):
        return "iso-dC/iso-dG"
    key = _END_PREFIX.sub("", key).lower().lstrip(" -")
    for rx, name in _CANON_RX:
        if rx.fullmatch(key):
            return name
    return ""


def _name_mods(end: str, text: str) -> list[str]:
    out = []
    for part in re.findall(rf"{_MOD}", text):
        part = part.strip(" -\t")
        if not part:
            continue
        c = canonical(part)
        shown = part[1:-1] if part[:1] == "(" and part[-1:] == ")" else part
        out.append(f"{end}' {c} ({shown})" if c and c != shown.lower() else f"{end}' {shown}")
    return out


def modifications(orig: str, start: int, end: int, seq: str = "") -> tuple[int, int, list[str]]:
    """Widen [start, end) in the original text over the 5' and 3' modifications touching
    it -> (start, end, a readable list of every modification found)."""
    mods = []
    if start > 0 and orig[start - 1] in "rm+" and orig[start] in BASE:
        start -= 1                        # the first base's own mark: rGrUrU...
    bases = orig[start:end]
    head = orig[max(0, start - 60):start].split("\n")[-1]
    m5 = MOD5.search(head) if head.strip(" -") and not head.endswith("\t") else None
    txt = re.sub(r"^5['′’]\s*-?\s*", "", m5.group(0)).strip(" -\t") if m5 else ""
    if txt:
        mods += _name_mods("5", txt)
        start -= len(head) - m5.start() - (len(m5.group(0)) - len(m5.group(0).lstrip()))
    tail = orig[end:end + 60].split("\n")[0]
    m3 = MOD3.match(tail) if tail.strip(" -") and not tail.startswith("\t") else None
    if m3 and m3.group(0).strip(" -"):
        mods += _name_mods("3", m3.group(0).strip(" -"))
        end += m3.end()
    if _ANY_MARK.search(bases):
        counts = {"r": len(_R_MARK.findall(bases)), "m": len(_M_MARK.findall(bases)),
                  "+": len(_L_MARK.findall(bases)),
                  "*": bases.count("*") + bases.count("∗") + bases.count("-s-")}
        mods += [f"{_TAG[k]} x{n}" for k, n in counts.items() if n]
        mods += [f"internal /{x}/" for x in IDT_INTERNAL.findall(bases)]
        mods += [f"internal {x}" for x in re.findall(r"\(d[UI]\)", bases)]
    else:
        counts = {"r": 0}
    if (re.search(r"3['′’]\s*[-–]?\s*\Z", orig[max(0, start - 6):start])
            and re.match(r"\s*[-–]?\s*5['′’]", orig[end:end + 6])):
        mods.append("written 3'->5' (reverse it to read 5'->3')")
    u = seq.count("U") + seq.count("u")
    if u:
        u -= len(re.findall(r"(?i)\(dU\)|/i?(?:deoxy|d)U/", bases))
    if u and not counts["r"]:
        mods.append("RNA (written with U)" if not re.search("[Tt]", seq) and u >= 2
                    else f"U bases x{u}")
    return start, end, mods


# ------------------------------------------------------------ the scan copy

class Scan:
    """A normalised copy of a text (marks removed, shorthand expanded, placeholders made
    bases or "…") with a map from every copy position back to a span of the original."""

    def __init__(self, orig: str, edits: list[tuple[int, int, str]]):
        self.orig = orig
        parts: list[str] = []
        self._scan: list[int] = []
        self._orig: list[int] = []
        self._olen: list[int] = []
        self._copy: list[bool] = []
        p = at = 0
        for a, b, r in sorted(edits) + [(len(orig), len(orig), "")]:
            for txt, o, ol, cp in ((orig[at:a], at, a - at, True), (r, a, b - a, False)):
                if txt:
                    self._scan.append(p)
                    self._orig.append(o)
                    self._olen.append(ol)
                    self._copy.append(cp)
                    parts.append(txt)
                    p += len(txt)
            at = b
        self.text = "".join(parts)
        self._nl: list[int] | None = None

    def span(self, a: int, b: int) -> tuple[int, int]:
        """The original span [start, end) that copy positions [a, b) came from."""
        k = bisect_right(self._scan, a) - 1
        s = self._orig[k] + (a - self._scan[k] if self._copy[k] else 0)
        k = bisect_right(self._scan, b - 1) - 1
        e = self._orig[k] + (b - self._scan[k] if self._copy[k] else self._olen[k])
        return s, e

    def line(self, o: int) -> int:
        """1-based line of original position o."""
        if self._nl is None:
            self._nl = [m.start() for m in re.finditer("\n", self.orig)]
        return bisect_left(self._nl, o) + 1


_BR = re.compile(r"[\[{<]([^\[\]{}<>\n]{1,40}(?:\n[^\[\]{}<>\n]{1,40})?)[\]}>]")
_BR_WORD = re.compile(r"(?i)barcode|\bbc\d*\b|\bcbc\b|\bubc\b|umi|index|\bidx|\bi[57]\b|"
                      r"\btag\b|\bwell\b|\bsample|round|random|insert|cdna|\bn\d|\bx\d|"
                      r"\bn-?mer\b|\d\s*-?\s*(?:bp|nt|bases?|mer)\b")
# (the patterns below start with a character, not a look-behind, so that re can skip
# ahead fast in a 50 MB table; what must not precede a match is checked in scan_copy)
_HOMO_PAREN = re.compile(r"\(d?([ACGTUN])\)[ ]?(\d{1,3})(?![0-9])")       # (T)30, (dT)30
_HOMO = re.compile(r"([ACGTUN])(?:\((\d{1,3})\)|(\d{1,3})(?![0-9]))")       # T(30), T30, N8
_PAREN_GAP = re.compile(r"\((?:d|p|poly)?[ACGTU]\)n?|\((?:d|p)?[ACGTU]\)[ ]?n\b|poly-?\(?[AT]\)?n\b")
_XRUN = re.compile(rf"XXX+(?:(?:\.{{2,}}|{GAP})XXX+)?(?![EFILOPQZa-qs-z])|JJJ+(?:(?:\.{{2,}}|{GAP})JJJ+)?"
                   rf"(?![EFILOPQZa-qs-z])")
_ELLIPSIS = re.compile(rf"\.\.\.(?=[{BASE}])|…(?=[{BASE}])")
_GLUED_DNA = re.compile(rf"[cgm](?:DNA|RNA)(?:_?[Ii]nsert)?")
_ISO = re.compile(r"i[CG]i[CG](?:i[CG])*")
_THIO_FAST = re.compile(r"[*∗]|-s-")
_IDT_FAST = re.compile(r"/(i?[A-Za-z][\w-]*)/")
_PARENS_BASE = re.compile(r"\(d([UI])\)")
_LINENO = re.compile(rf"\n[ \t\f]*\n[ \t\f]*\d{{1,4}}[ \t]*\n[ \t\f]*\n[ \t\f]*"
                     rf"(?=[{BASE}{LOW}]{{4}})")
_WRAP = re.compile(rf"[{BASE}{LOW}][ ]*\n ?(?=[{BASE}{LOW}])")
_INDENT = re.compile(r"\n(?=[ \t]{2,}\S)")
_SEQLINE = re.compile(rf"(?m)^[ \t]*([{BASE}{LOW}]{{8,}})[ \t]*\r?$")


def _bracket(m: re.Match) -> str | None:
    """What a [placeholder] stands for: its bases, N x its length, or "…"; None if it is
    not a placeholder ([1-32], [5Phos], a citation)."""
    body = m.group(1).replace("\n", " ").strip()
    if re.fullmatch(f"[{BASE}]{{3,}}", body) and not canonical(body):
        return body                                 # [NNNNNN]
    h = re.fullmatch(r"(?:poly-?)?d?\(?([ACGTUN])\)?(\d{1,3})", body)
    if h and 3 <= int(h.group(2)) <= 200:
        return h.group(1) * int(h.group(2))
    if not _BR_WORD.search(body):
        return None
    if re.search(r"\bor\b|\d\s*[-–]\s*\d|/", body):
        return GAP                                  # "9 bp or 10 bp", "8-10 nt"
    n = {int(x) for x in re.findall(r"(\d+)\s*-?\s*(?:bp|nt|bases?|mer)\b", body, re.I)}
    n |= {int(x) for x in re.findall(r"(?i)(?:umi|bc|cbc|ubc|index|idx|barcode|\bn|\bx)[ _-]?(\d+)\b", body)}
    n = {x for x in n if 1 <= x <= 60}
    return "N" * n.pop() if len(n) == 1 else GAP


def scan_copy(orig: str) -> Scan:
    """The normalised copy the scanner and --find read; see the module docstring."""
    edits: list[tuple[int, int, str]] = []
    taken = bytearray(len(orig) + 1)       # 1 where an earlier (higher-priority) edit sits

    def add(a: int, b: int, r: str) -> bool:
        if a >= b or taken.find(1, a, b) >= 0:
            return False
        taken[a:b] = b"\1" * (b - a)
        edits.append((a, b, r))
        return True

    for m in PER_BASE.finditer(orig):
        add(m.start(), m.end(), re.sub(r"[rm+*∗\s]", "", m.group(0)))
    def prev(a: int) -> str:
        return orig[a - 1] if a else ""

    for m in _IDT_FAST.finditer(orig):
        if not (prev(m.start()) in BASE + LOW + "]" and prev(m.start())
                and orig[m.end():m.end() + 1] in tuple(BASE + LOW + "[")):
            continue
        name = m.group(1).lower()
        add(m.start(), m.end(), "U" if re.search(r"deoxyu|^i?du$", name) else
            "N" if re.search(r"deoxyi|^i?di$|inosine", name) else "")
    for m in _PARENS_BASE.finditer(orig):
        add(m.start(), m.end(), "U" if m.group(1) == "U" else "N")
    for m in _BR.finditer(orig):
        r = _bracket(m)
        if r is not None:
            add(m.start(), m.end(), r)
    for m in _XRUN.finditer(orig):
        g = m.group(0)
        if re.match(r"[EFILOPQZa-z]", prev(m.start())):
            continue
        add(m.start(), m.end(), GAP if "." in g or GAP in g else "N" * len(g))
    for m in _ELLIPSIS.finditer(orig):
        if prev(m.start()) and prev(m.start()) in BASE:
            add(m.start(), m.end(), GAP)
    for m in _PAREN_GAP.finditer(orig):
        a, b = m.start(), m.end()
        if (re.search(rf"[{BASE}{LOW}{GAP}\]][-]?\Z", orig[max(0, a - 2):a])
                or re.match(rf"-?[{BASE}{LOW}\[]", orig[b:b + 2])):
            add(a, b, GAP)
    # homopolymer shorthand (T30, dT30, (T)30, T(30), N8) only next to bases or another
    # placeholder, so "T30 oligo" and "A549 cells" in prose stay prose
    cand = []
    for m in list(_HOMO_PAREN.finditer(orig)) + list(_HOMO.finditer(orig)):
        paren = m.re is _HOMO_PAREN
        base, n = (m.group(1), m.group(2)) if paren else (m.group(1), m.group(2) or m.group(3))
        bare = not paren and m.group(3) is not None
        a0 = m.start() - (not paren and prev(m.start()) == "d")     # dT30, dT(30)
        p = prev(a0)
        if p and (p.isdigit() or p.islower() or p in "EFIJLOPQXZ" if bare else p.isalnum()):
            continue
        least = 10 if bare and base != "N" else 3      # bare "T7", "U6" are names
        if bare and base == "U":
            continue
        if least <= int(n) <= 200:
            cand.append((a0, m.end(), base * int(n)))
    done = set()
    for _ in range(4):
        changed = False
        for a, b, r in cand:
            if (a, b) in done:
                continue
            left = orig[max(0, a - 5):a]
            right = orig[b:b + 5]
            ok = (re.search(rf"[{BASE}{LOW}]{{4}}-?\Z", left) or re.match(rf"-?[{BASE}]{{4}}", right)
                  or (r[0] in "AT" and re.match(r"[VN]{1,2}(?![A-Za-z])", right))
                  or (a and taken[a - 1]) or (a > 1 and orig[a - 1] == "-" and taken[a - 2])
                  or taken[b] or (orig[b:b + 1] == "-" and taken[b + 1]))
            if ok and add(a, b, r):
                done.add((a, b))
                changed = True
        if not changed:
            break
    for m in _ISO.finditer(orig):
        if not prev(m.start()).isalpha():
            add(m.start(), m.end(), "#" * (m.end() - m.start()))
    for m in _GLUED_DNA.finditer(orig):
        a, b = m.start(), m.end()
        if not (prev(a) and prev(a) in BASE + LOW + "-"):
            continue
        between = (re.search(rf"[{BASE}]-?\Z", orig[max(0, a - 2):a])
                   and re.match(rf"-?[{BASE}]{{4}}", orig[b:b + 5]))
        add(a, b, GAP if between else "#" * (b - a))      # "...GGG-cDNA_Insert-AGATC..."
    for m in _THIO_FAST.finditer(orig):
        a, b = m.start(), m.end()
        if m.group(0) != "-s-":         # "G* A*C", "G*A\n*T": the mark may sit in a wrap
            while a and orig[a - 1] in " \t":
                a -= 1
            a -= orig[a - 1:a] == "\n"
            while a and orig[a - 1] in " \t":
                a -= 1
            stop = re.match(r"[ \t]*\n?[ \t]*", orig[b:b + 40])
            b += stop.end()
        if prev(a) and prev(a) in BASE + LOW and orig[b:b + 1] and orig[b] in BASE + LOW:
            add(a, b, "")
    for m in _LINENO.finditer(orig):
        if not (prev(m.start()) and prev(m.start()) in BASE + LOW):
            continue
        line_before = orig[orig.rfind("\n", 0, m.start()) + 1:m.start()]
        if len(line_before) >= 30:
            add(m.start(), m.end(), "\n")
    # "...CCG*A\n*T*C*T\nCAAGCAG...": a phosphorothioate 3' end closes the oligo
    for m in re.finditer(rf"(?<=[*∗][{BASE}{LOW}])[ \t]*\n(?![ \t]*[*∗])", orig):
        add(m.end() - 1, m.end(), HARD_NL)
    # a spreadsheet or Word table row never wraps: no join out of or into a tabbed line
    for m in _WRAP.finditer(orig):
        nl = m.start() + m.group(0).index("\n")
        ls = orig.rfind("\n", 0, nl) + 1
        le = orig.find("\n", nl + 1)
        if "\t" in orig[ls:nl] or "\t" in orig[nl + 1:le if le >= 0 else len(orig)]:
            add(nl, nl + 1, HARD_NL)
    for m in _INDENT.finditer(orig):
        add(m.start(), m.start() + 1, HARD_NL)
    # 3+ consecutive lines that are each one sequence of the same length: a list
    prev = None
    block: list[re.Match] = []
    for m in list(_SEQLINE.finditer(orig)) + [None]:
        if m is not None and prev is not None and m.start() == prev.end() + 1 \
                and len(m.group(1)) == len(block[0].group(1)):
            block.append(m)
        else:
            if len(block) >= 3:
                for x in block[:-1]:
                    add(x.end(), x.end() + 1, HARD_NL)
            block = [m] if m is not None else []
        prev = m
    for m in re.finditer(r"\r(?=\n)", orig):
        add(m.start(), m.end(), "")
    return Scan(orig, edits)


def unmodify(text: str) -> tuple[str, list[int]]:
    """-> (the scan copy, the original position of each of its characters). Kept for
    callers of the old interface; the scanner itself uses scan_copy()."""
    sc = scan_copy(text)
    return sc.text, [sc.span(i, i + 1)[0] for i in range(len(sc.text))]


# --------------------------------------------------------------------- scanning

_L = f"[{BASE}{LOW}{GAP}]"
_UP = f"[{BASE}{GAP}]"
_SEP = r"(?:[ \t]*\n ?| |-)"
# The first base is matched before the look-behinds that vet it, so re can skip ahead
# fast. A run may not continue a word or number; a lower-case base may not follow any
# letter ("Biotin" + "ACAC..." is not "nACAC..."), an upper-case one may follow a
# modification glued to it (rAppTGG..., /5Phos/ACG...).
_FIRST = (rf"[{BASE}{LOW}{GAP}](?<![A-Z0-9{GAP}].)"
          rf"(?:(?<=[{BASE}{GAP}])(?<![^A-Za-z][{LOW}].)|(?<![a-z].))")
_FIRST_UP = rf"[{BASE}](?<![A-Z0-9{GAP}].)(?<![^A-Za-z][{LOW}].)"
_END = r"(?=[ \t]*(?:[\t\n\v]|$)|[rm+*∗/]|\s*[-–]\s*[35]['′’])"
RUN = re.compile(rf"(?:{_FIRST_UP}[{BASE}]{{0,2}} (?:[{BASE}]{{3}} ){{3,}}{_L}+|{_FIRST}{_L}{{3,}})"
                 rf"(?:{_SEP}(?:{_L}{{4,}}|{_UP}{{1,3}}(?={_SEP}{_L}{{4}}|{_END})))*"
                 rf"(?![A-Z0-9{LOW}])")

_NAME = re.compile(r"[A-Za-z0-9](?:[\w.+/'′–-]|\[[\w.,–-]{1,12}\])*")
_NOT_NAME = re.compile(r"^(?:NH2|rApp|App|biotin|Bio|Btn|Phos|5Phos|5['′]|3['′]|5|3|Fig|Table|"
                       r"Supplementary|DNA|RNA|cDNA|PCR|ml|µl|μl|nM|μM|µM|UMI|UMIs|barcodes?|"
                       r"\d+-?(?:bp|nt|mer)|bp|nt|[NX]\d+|(?:i[CG])+i?|sequence|Sequence|Oligo|"
                       r"oligo|primer|Primer|X+|J+|([A-Z])\1{3,})$", re.I)


_TRAILING_MOD = re.compile(r"\[[^\]]*\]$|/[\w-]+/\s*$")
_OLIGO_TOKEN = re.compile(f"[{BASE}acgtun*+X.…]{{8,}}")
_BARCODE_TOKEN = re.compile(r"[ACGTN]{4,7}")


def label_before(before: str) -> str:
    """The oligo's name: the nearest token before it that looks like a name -- has a
    capital or a digit ("RA3 (rApp...", "P5 | AATG...", "Oligo-dT30VN: 5'-...") --
    skipping the modification written in front of the bases and a barcode column."""
    before = _TRAILING_MOD.sub("", before.rstrip())
    for tok in reversed(_NAME.findall(before[-120:])[-5:]):
        tok = tok.strip("-'′.:")
        if _OLIGO_TOKEN.fullmatch(tok):
            return ""                     # an earlier oligo: its name is not this one's
        if _BARCODE_TOKEN.fullmatch(tok) or _NOT_NAME.match(tok) or canonical(tok):
            continue                      # a barcode cell, a unit, a modification
        if len(tok) >= 2 and (any(c.isupper() for c in tok) or any(c.isdigit() for c in tok)):
            return tok
    return ""


def _row_label(row_before: str) -> str:
    """The name of a table row: of the cells before the oligo's, the nearest one that is
    a name; an identifier (no spaces: TN5_A_ME, v3_gr1) is preferred to a description
    ("Indexed Tn5 Transposome assembly"), which is kept only if there is no identifier."""
    described = ""
    for c in reversed(row_before.split("\t")[:-1]):
        c = c.strip()
        if (not c or len(c) > 60 or not re.search("[A-Za-z]", c) or re.fullmatch(r"[A-P]\d{1,2}", c)
                or re.fullmatch(f"[{BASE}acgtun]+", c)):
            continue
        if " " not in c:
            return c
        described = described or c
    return described


def context_at(text: str, start: int, end: int, width: int = 400) -> str:
    """For a table row (a tab on its line): that row, cells joined by ' | '. Otherwise the
    sentence around [start, end): back to the previous '. ' or blank line, forward to the
    next one, clipped to `width` either side."""
    ls = text.rfind("\n", 0, start) + 1
    le = text.find("\n", end)
    le = len(text) if le < 0 else le
    if "\t" in text[ls:le] and le - ls < 4000:
        return re.sub(r"[ \t]*\t[ \t]*", " | ", text[ls:le]).strip(" |")
    return sentence_at(text, start, end, width)


def sentence_at(text: str, start: int, end: int, width: int = 400) -> str:
    lo = max(0, start - width)
    left = text[lo:start]
    cut = max(left.rfind(". "), left.rfind("\n\n"), left.rfind(".\n"))
    left = left[cut + 2:] if cut >= 0 else left
    right = text[end:end + width]
    m = re.search(r"\.\s|\n\s*\n", right)
    right = right[:m.start() + 1] if m else right
    return re.sub(r"\s+", " ", left + text[start:end] + right).strip()


_SPACE_DASH = re.compile(r"[\s\-]")
_NRUN = re.compile(r"N{3,}")
_WS = re.compile(r"\s+")
_BR_SPACE = re.compile(r"\[[^\[\]]{1,60}\]|\{[^{}]{1,60}\}|\([^()]{1,30}\)")
_LINENO_TEXT = re.compile(r"\n[ \t\f]*\n[ \t\f]*\d{1,4}[ \t]*\n")


def scan(path: Path, orig: str, known: dict[str, str], min_len: int = MIN_LEN,
         width: int = 400, sc: Scan | None = None) -> list[Hit]:
    sc = sc or scan_copy(orig)
    text = sc.text
    hits: list[Hit] = []
    seen: dict[str, Hit] = {}
    for m in RUN.finditer(text):
        shown = _SPACE_DASH.sub("", m.group(0))
        seq = shown.upper().replace("U", "T")
        body = seq.replace(GAP, "")
        concrete = sum(body.count(b) for b in "ACGT")
        rest = _NRUN.sub("", body) if "NNN" in body else body
        if (len(body) < min_len or concrete < MIN_CONCRETE
                or concrete / max(1, len(rest)) < MIN_ACGT):
            continue
        s0, e0 = sc.span(m.start(), m.end())
        s, e, mods = modifications(orig, s0, e0, shown)
        ls = orig.rfind("\n", 0, s) + 1
        le = orig.find("\n", e)
        cell = (not orig[ls:s].rsplit("\t", 1)[-1].strip(" ")
                and not orig[e:le if le >= 0 else len(orig)].split("\t", 1)[0].strip(" "))
        if (shown != seq and sum(c.isupper() for c in shown) < MIN_UPPER and not cell
                and not mods):
            continue
        before = orig[max(0, s - 80):s].split("\n")[-1]
        label = _row_label(orig[ls:s]) if cell and "\t" in orig[ls:s] else label_before(before)
        if not label and not before.strip():          # the name on the line above
            above = orig[max(0, s - 200):s].rstrip().split("\n")[-1].strip()
            above = re.sub(r"^\d{1,2}[.)]\s+|[\s:]+$", "", above)
            if (2 <= len(above) <= 50 and label_before(above) and not re.search("[ACGTU]{8}", above)
                    and not re.fullmatch(r"[\d.,\s]+\S{0,3}|.*\b(?:Table|Figure|Fig)\b.*", above)):
                label = above if "\t" not in above else label_before(above)
        line = sc.line(s)
        written = orig[s:e]
        if "\n" in written or " " in written or "\t" in written:
            written = _WS.sub("", _BR_SPACE.sub(lambda b: _WS.sub("\x00", b.group(0)),
                                                _LINENO_TEXT.sub("\n", written)))
            written = written.replace("\x00", " ")
        first = seen.get(written)
        if first:                         # scg pages draw one oligo in many panels
            first.also.append(line)
            continue
        seen[written] = h = Hit(path, line, shown, context_at(orig, s, e, width),
                                label, contained(seq, known), written, mods,
                                table_cell=cell)
        hits.append(h)
    return hits


# --------------------------------------------------------------- families

def _stem(label: str) -> str:
    """'P7-1' -> 'P7', 'sc_ligation_12' -> 'sc_ligation', 'MDRT001' -> 'MDRT',
    'HYi7_1_CGCTCAGTTC' -> 'HYi7': the name of a table's rows without the row number."""
    s = re.sub(r"[\W_]+[ACGTN]{4,}$", "", label)
    return re.sub(r"[\W_]*(?:[A-Pa-p]?\d+)$", "", s)


def collapse(hits: list[Hit], min_size: int = 4) -> list[Hit]:
    """Fold runs of >= min_size consecutive hits that are one oligo with a varying
    barcode or index into the first of them, with .family and .windows set. Rows belong
    together when they have the same end modifications and either the same name stem
    (P7-1, P7-2, ...) and a common 12-nt stretch, or one length and differences at only a
    minority of positions."""
    def ends(h: Hit) -> list[str]:
        return [x for x in h.mods if x[:2] in ("5'", "3'")]

    def kmers(x: str) -> set[str]:
        return {x[k:k + 12] for k in range(len(x) - 11)}

    out: list[Hit] = []
    i = 0
    while i < len(hits):
        h = hits[i]
        ref = h.seq.upper()
        n = len(ref)
        stem = _stem(h.label)
        ref12 = kmers(ref)
        vary: set[int] = set()
        j = i + 1
        while j < len(hits):
            g = hits[j].seq.upper()
            if GAP in g or ends(hits[j]) != ends(h):
                break
            if len(g) == n:
                v = vary | {k for k in range(n) if g[k] != ref[k]}
                if len(v) <= max(12, n // 2) and n - len(v) >= 12:
                    vary = v
                    j += 1
                    continue
            if stem and _stem(hits[j].label) == stem and ref12 & kmers(g):
                j += 1
                continue
            break
        fam = hits[i:j]
        if len(fam) >= min_size:
            h.family = fam
            seqs = [x.seq.upper() for x in fam]
            if len({len(x) for x in seqs}) == 1:
                pos = sorted({k for x in seqs for k in range(n) if x[k] != ref[k]})
            else:
                p = min(len(_common_prefix(ref, x)) for x in seqs)
                q = min(len(_common_prefix(ref[::-1], x[::-1])) for x in seqs)
                pos = []
                h.windows = [(p + 1, max(p + 1, n - q))]
            runs = [[pos[0], pos[0]]] if pos else []
            for k in pos[1:]:
                if k - runs[-1][1] <= 3:
                    runs[-1][1] = k
                else:
                    runs.append([k, k])
            if pos:
                h.windows = [(a + 1, b + 1) for a, b in runs]
            i = j
        else:
            i += 1
        out.append(h)
    return out


def _common_prefix(a: str, b: str) -> str:
    k = 0
    while k < min(len(a), len(b)) and a[k] == b[k]:
        k += 1
    return a[:k]


# ------------------------------------------------------------------------- find

_FIND_GAP = r"(?:[ \t\n\-*∗]|-s-|/[\w-]{1,20}/|\+?[rm+](?=[A-Z]))*"


_MODLIKE = re.compile(rf"(?:{_MOD})+")


_TRIM_GAP = re.compile(rf"(.*[{BASE}{LOW}{GAP}])(?:{_FIND_GAP})", re.S)


def _letters(q: str) -> str:
    """The bases of a query as pasted: ends, /mods/, marks, shorthand and placeholders
    handled like the text; raises ValueError on letters that are not bases."""
    q = re.sub(rf"(?:(?<=[{BASE}*∗ '′’-])|^)(?:\+[rm]|[rm+])(?=[{BASE}])", "", q.strip())  # rG, +G, +rG
    sc = scan_copy(" " + q + " ")
    runs = list(RUN.finditer(sc.text))
    if not runs:
        raise ValueError(f"no bases in {q!r}")
    t = sc.text[runs[0].start():runs[-1].end()]   # 5'/3' ends and their modifications dropped
    for end, rx in ((sc.text[:runs[0].start()], MOD5), (sc.text[runs[-1].end():], MOD3)):
        rest = re.sub(r"^\s*(?:5|3)['′’]?\s*-?|-?\s*(?:5|3)['′’]?\s*$", "", end).strip(" -")
        if rest and not (rx.fullmatch(rest) or rx.fullmatch(rest + "-")
                         or all(canonical(x) or _MODLIKE.fullmatch(x) for x in re.split(r"[\s-]+", rest) if x)):
            raise ValueError(f"not bases or a known modification: {rest!r} in {q!r}")
    t = re.sub(r"[\s\-*∗#]|/[\w-]+/", "", t)
    t = re.sub(rf"[rm+](?=[{BASE}])", "", t)
    bad = sorted(set(re.sub(f"[{BASE}{LOW}{GAP}]", "", t)))
    if bad:
        raise ValueError(f"not bases: {''.join(bad)!r} in {t!r}")
    return t.upper().replace("U", "T")


def _pattern(q: str, strict: bool) -> str:
    parts = []
    for run in re.finditer(rf"N+|{GAP}|.", q):
        r = run.group(0)
        if r == GAP:
            parts.append(rf"(?:{GAP}|(?:[{BASE}{LOW}]{_FIND_GAP}){{1,60}}?)")
            continue
        c = r[0]
        if strict:
            cls = "[TtUu]" if c == "T" else f"[{c}{c.lower()}]"
        else:
            ok = [t for t in _IUPAC if set(_IUPAC[t]) & set(_IUPAC[c])]
            cls = "[" + "".join(ok) + "".join(t.lower() for t in ok if t in "ACGTUN") + "]"
        unit = f"{cls}{_FIND_GAP}"
        if c == "N":
            parts.append(rf"(?:(?:{unit}){{{len(r)}}}|{GAP}{_FIND_GAP})")
        else:
            parts.append(unit * len(r))
    return "".join(parts)


def find(text: str, query: str, strict: bool = False, sc: Scan | None = None,
         detail: bool = False) -> list:
    """Line numbers where `query` occurs on either strand (or, with detail=True, tuples
    (line, strand, 'exact'|'iupac', as written)); see the module docstring."""
    q = _letters(query)
    if not q:
        return []
    sc = sc or scan_copy(text)
    seen = set()
    out = []
    orients = [("+", q), ("rc", revcomp(q)), ("rev", q[::-1]), ("comp", complement(q))]
    done = set()
    for strand, x in orients:
        if x in done:
            continue
        done.add(x)
        rx = re.compile(_pattern(x, strict))
        for m in rx.finditer(sc.text):
            mtext = _TRIM_GAP.fullmatch(m.group(0)).group(1)
            s, e = sc.span(m.start(), m.start() + len(mtext))
            line = sc.line(s)
            if (line, s) in seen:
                continue
            seen.add((line, s))
            got = re.sub(rf"[rm+](?=[A-Z])|/[\w-]+/|-s-|[^{BASE}{LOW}{GAP}]", "", mtext)
            got = got.upper().replace("U", "T")
            exact = got == x
            if not exact and not strict and len(got) == len(x):
                # mostly placeholders in the text (a drawn NNNN...NNNN) is not a location
                fixed = [i for i, c in enumerate(x) if c in "ACGT"]
                if sum(got[i] == x[i] for i in fixed) < 0.6 * len(fixed):
                    continue
            written = re.sub(r"[ \t]*\n[ \t]*", "", sc.orig[s:e])
            out.append((line, strand, "exact" if exact else "iupac", written))
    out.sort(key=lambda t: (t[2] != "exact", t[0]))
    return out if detail else sorted({t[0] for t in out})


# ------------------------------------------------------------------------ main

def _read_all(paths: list[str]):
    """(path, text) for every readable file, skipping one whose content was already read."""
    have: dict[str, Path] = {}
    for path in files(paths):
        size = path.stat().st_size
        if size > 20e6:
            print(f"reading {path} ({size / 1e6:.0f} MB) ...", file=sys.stderr)
        text = read(path)
        if text is None:
            continue
        digest = hashlib.sha1(text.encode("utf-8", "replace")).hexdigest()
        if digest in have:
            print(f"skip {path}: same content as {have[digest]}", file=sys.stderr)
            continue
        have[digest] = path
        yield path, text


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                 epilog=__doc__.split("\n\n", 1)[1],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="+", help="files or directories")
    ap.add_argument("--min-len", type=int, default=MIN_LEN, help="shortest sequence kept (nt)")
    ap.add_argument("--known-only", action="store_true",
                    help="only hits containing a sequence already in lib/")
    ap.add_argument("--keyword-only", action="store_true",
                    help="only hits whose sentence mentions primer/oligo/adapter/...")
    ap.add_argument("--find", metavar="SEQ", action="append",
                    help="instead: where does SEQ occur (either strand); repeatable")
    ap.add_argument("--strict", action="store_true",
                    help="--find: letters must be equal (N matches only N), not IUPAC classes")
    ap.add_argument("--tsv", action="store_true", help="tab-separated output")
    ap.add_argument("--max-hits", type=int, default=200,
                    help="per file: hits (after collapsing families) or --find locations "
                         "shown; the rest are counted. 0 = no limit")
    ap.add_argument("--all", action="store_true",
                    help="list every hit: no family collapsing and no --max-hits cut")
    ap.add_argument("--context", type=int, default=300, metavar="CHARS",
                    help="context characters either side of a hit outside a table (300)")
    a = ap.parse_args(argv)
    limit = 0 if a.all else a.max_hits

    if a.find:
        try:
            queries = [(f, _letters(f)) for f in a.find]
        except ValueError as e:
            print(f"--find: {e}", file=sys.stderr)
            return 2
        if a.tsv:
            print("query\tfile\tline\tstrand\tmatch\twritten")
        totals = {f: [0, 0] for f, _ in queries}
        for path, text in _read_all(a.paths):
            sc = scan_copy(text)
            for f, _ in queries:
                locs = find(text, f, a.strict, sc, detail=True)
                for i, (ln, strand, kind, written) in enumerate(locs):
                    totals[f][0] += 1
                    totals[f][1] += kind == "exact"
                    if limit and i >= limit:
                        rest = locs[i:]
                        if not a.tsv:
                            print(f"{path}: ... {len(rest)} more location(s), lines "
                                  f"{min(t[0] for t in rest)}-{max(t[0] for t in rest)}; --max-hits 0 shows all")
                        totals[f][0] += len(rest) - 1
                        totals[f][1] += sum(k == "exact" for _, _, k, _ in rest[1:])
                        break
                    w = written if len(written) <= 120 else written[:117] + "..."
                    if a.tsv:
                        print("\t".join([f, str(path), str(ln), strand, kind, written]))
                    else:
                        tag = f"  [{f}]" if len(queries) > 1 else ""
                        print(f"{path}:{ln}  {strand:<4} {kind:<5}  {w}{tag}")
        for f, (n, ex) in totals.items():
            print(f"{n} location(s) of {f}: {ex} exact, {n - ex} only through IUPAC letters",
                  file=sys.stderr)
        return 0

    known = known_sequences()
    # Read everything first, so the files can be ordered by their most oligo-like hit:
    # the primer table of a supplement comes before its motif-enrichment table.
    per_file = []
    for path, text in _read_all(a.paths):
        hits = [h for h in scan(path, text, known, a.min_len, a.context)
                if (h.known or not a.known_only) and (h.keyword or not a.keyword_only)]
        if not a.all:
            hits = collapse(hits)
        dropped = 0
        if limit and len(hits) > limit:   # keep the most oligo-like, in file order
            keep = {id(h) for h in sorted(hits, key=lambda h: -h.score)[:limit]}
            dropped = len(hits) - limit
            hits = [h for h in hits if id(h) in keep]
        if hits:
            per_file.append((max(h.score for h in hits), path, hits, dropped))
    per_file.sort(key=lambda t: -t[0])

    if a.tsv:
        print("file\tline\tlength\tscore\tlabel\tsequence\twritten\tmodifications\tknown\t"
              "keyword\talso_at\tsentence\tfamily_size\tfamily_lines\tvariable_nt")
    for _, path, hits, dropped in per_file:
        for h in hits:
            fam = (f"{len(h.family)}", f"{h.family[0].line}-{h.family[-1].line}",
                   ", ".join(f"{a}-{b}" for a, b in h.windows)) if h.family else ("", "", "")
            if a.tsv:
                print("\t".join([str(h.path), str(h.line), h.length_text.split()[0],
                                 str(h.score), h.label, h.seq, h.written, "; ".join(h.mods),
                                 "; ".join(h.known), "yes" if h.keyword else "",
                                 ",".join(map(str, h.also)), h.sentence, *fam]))
                continue
            lab = f"{h.label}: " if h.label else ""
            out = [f"{h.path}:{h.line}  {lab}{h.seq} ({h.length_text}, score {h.score})"]
            if h.written != h.seq:
                out.append(f"    as written: {h.written}")
            if h.mods:
                out.append(f"    modifications: {'; '.join(h.mods)}")
            if h.known:
                out.append(f"    contains: {', '.join(h.known)}")
            if h.family:
                lens = sorted({x.length for x in h.family})
                ln = f" ({lens[0]}-{lens[-1]} nt long)" if len(lens) > 1 else ""
                out.append(f"    family: {len(h.family)} variants, lines {fam[1]}, variable at "
                           f"nt {fam[2]}{ln} (this is the first; --all lists every row)")
            if h.also:
                out.append(f"    also at lines: {', '.join(map(str, h.also))}")
            out.append(f"    context: {h.sentence}")
            print("\n".join(out) + "\n")
        if dropped and not a.tsv:
            print(f"{path}: ... {dropped} lower-scoring sequences not shown -- often a motif "
                  f"table or a sequence library; --max-hits 0 or --all shows them\n")
    print(f"{sum(len(t[2]) for t in per_file)} sequence(s) in {len(per_file)} file(s)",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
