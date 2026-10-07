#!/usr/bin/env python3
"""Turn source documents into plain text: each FILE gets a sibling FILE.txt.

Reading a paper means grepping it, quoting line numbers and handing it to an LLM, and all
of that is simpler on text than on PDF, Word, Excel or HTML. So every source is converted
once, next to the original, and everything downstream reads the .txt:

    python3 tools/doctext.py _data/sources/<slug>/          # a directory: every document in it
    python3 tools/doctext.py paper.pdf supp.docx            # -> paper.pdf.txt, supp.docx.txt
    python3 tools/doctext.py paper.pdf --stdout             # print instead of writing

Formats: .pdf (poppler's pdftotext, reading order), .docx (paragraphs; table cells
tab-separated), .xlsx (every sheet, tab-separated, under '## sheet <name>'), .html/.htm
(tags stripped, block elements become line breaks), .xml (JATS full text from PMC /
Europe PMC / eLife: tags stripped, paragraphs/titles/table rows on their own lines, cells
tab-separated), and -- when LibreOffice (`soffice`) is installed -- the legacy binary
.xls/.doc/.ppt (converted to xlsx/docx/pdf first, then read as those). A PDF with no text
layer (a scanned or "Print to PDF" image) is reported on stderr. A .txt is rewritten only
when its source is newer. tools/get_sources.py runs this on everything it downloads.
"""
from __future__ import annotations

import argparse
import html
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from scrape_primers import HTML_EXT, OFFICE_EXT  # noqa: E402
from scrape_primers import read as _read  # noqa: E402

XML_EXT = {".xml"}
LEGACY_EXT = {".xls", ".doc", ".ppt"}
SOFFICE = shutil.which("soffice") or shutil.which("libreoffice")
CONVERT = {".pdf"} | OFFICE_EXT | HTML_EXT | XML_EXT | (LEGACY_EXT if SOFFICE else set())

_BLOCK = ("p|title|article-title|subtitle|label|caption|tr|list-item|ref|abstract|sec|"
          "def-item|term|fn|table-wrap|fig|disp-quote|kwd|contrib|aff|mixed-citation|"
          "element-citation|statement|boxed-text|supplementary-material|disp-formula")

_INLINE = ("italic|bold|sup|sub|sc|underline|monospace|roman|sans-serif|overline|strike|"
           "styled-content|named-content|xref|ext-link|uri|inline-formula|email|abbrev")


def jats_text(raw: str) -> str:
    """JATS (or any article XML) -> text: one block element per line, table cells tab-
    separated, entities decoded. Front matter is kept (title, authors, abstract)."""
    raw = re.sub(r"(?s)<\?xml.*?\?>|<!DOCTYPE[^>]*>|<!--.*?-->", "", raw)
    raw = re.sub(r"(?s)<(mml:math|tex-math)\b.*?</\1>", " ", raw)
    raw = re.sub(r"(?i)</(?:td|th)>", "\t", raw)
    raw = re.sub(rf"(?i)</(?:{_BLOCK})>|<(?:{_BLOCK})\b[^>]*/>|<break\s*/>", "\n", raw)
    raw = re.sub(rf"(?i)<(?:{_BLOCK}|title-group|body|back|front|table)\b[^>]*>", "\n", raw)
    raw = re.sub(rf"</?(?:{_INLINE})\b[^>]*>", "", raw)     # markup inside a word/sequence
    raw = re.sub(r"<[^>]+>", " ", raw)                        # field boundaries (surname|given-names)
    text = html.unescape(raw)
    lines = [re.sub(r" *\t *", "\t", re.sub(r"[ \r]+", " ", ln)).strip(" ") for ln in text.split("\n")]
    out, blank = [], False
    for ln in lines:
        if ln.strip("\t "):
            out.append(ln.rstrip("\t "))
            blank = False
        elif not blank and out:
            out.append("")
            blank = True
    return "\n".join(out).strip() + "\n"


def _legacy(path: Path) -> str | None:
    """.xls/.doc/.ppt via LibreOffice -> xlsx/docx/pdf in a scratch dir, then read that."""
    if not SOFFICE:
        print(f"skip {path}: legacy {path.suffix} needs LibreOffice (soffice)", file=sys.stderr)
        return None
    target = {".xls": "xlsx", ".doc": "docx", ".ppt": "pdf"}[path.suffix.lower()]
    with tempfile.TemporaryDirectory() as tmp:
        profile = Path(tmp) / "profile"      # private profile: parallel runs do not collide
        done = subprocess.run([SOFFICE, f"-env:UserInstallation=file://{profile}", "--headless",
                               "--convert-to", target, "--outdir", tmp, str(path.resolve())],
                              capture_output=True, timeout=300, check=False)
        out = Path(tmp) / (path.stem + "." + target)
        if done.returncode or not out.exists():
            print(f"skip {path}: LibreOffice could not convert it", file=sys.stderr)
            return None
        return _read(out)


def read(path: Path) -> str | None:
    """Text of a document (None if unreadable). Extends scrape_primers.read with XML and
    the legacy Office formats."""
    ext = path.suffix.lower()
    if ext in XML_EXT:
        return jats_text(path.read_text(encoding="utf-8", errors="replace"))
    if ext in LEGACY_EXT:
        return _legacy(path)
    text = _read(path)
    if text is not None and ext == ".pdf" and len(text.strip()) < 100:
        print(f"warning {path}: no text layer (image-only PDF?) -- needs OCR or reading the "
              f"rendered pages", file=sys.stderr)
    return text


def text_path(src: Path) -> Path:
    return src.with_name(src.name + ".txt")


def convert(src: Path, force: bool = False) -> Path | None:
    """Write src.txt next to src (if stale) and return its path; None if unreadable."""
    out = text_path(src)
    if out.exists() and not force and out.stat().st_mtime >= src.stat().st_mtime:
        return out
    text = read(src)
    if text is None:
        return None
    out.write_text(text, encoding="utf-8")
    return out


def documents(paths: list[str]) -> list[Path]:
    found = []
    for a in paths:
        p = Path(a)
        if p.is_dir():
            found += sorted(f for f in p.rglob("*") if f.is_file() and f.suffix.lower() in CONVERT
                            and "_invalid" not in f.parts and ".git" not in f.parts)
        elif p.exists():
            found.append(p)
        else:
            print(f"skip {a}: no such file", file=sys.stderr)
    return found


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("paths", nargs="+", help="files or directories")
    ap.add_argument("--stdout", action="store_true", help="print the text, write nothing")
    ap.add_argument("--force", action="store_true", help="rewrite even if up to date")
    a = ap.parse_args(argv)
    for src in documents(a.paths):
        if a.stdout:
            text = read(src)
            if text is not None:
                sys.stdout.write(text)
            continue
        out = convert(src, a.force)
        if out:
            print(f"{out}  ({out.stat().st_size:,} bytes)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
