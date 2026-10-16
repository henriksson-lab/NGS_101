#!/usr/bin/env python3
"""Turn source documents into plain text: each FILE gets a sibling FILE.txt.

Reading a paper means grepping it, quoting line numbers and handing it to an LLM, and all
of that is simpler on text than on PDF, Word, Excel or HTML. So every source is converted
once, next to the original, and everything downstream reads the .txt:

    python3 tools/doctext.py /husky/henriksson/NGS101_papers/sources/<slug>/
    python3 tools/doctext.py paper.pdf supp.docx            # -> paper.pdf.txt, supp.docx.txt
    python3 tools/doctext.py paper.pdf --stdout             # print instead of writing

Formats: .eps (via Ghostscript when installed, then the ordinary PDF reader),
.pdf (poppler's pdftotext, reading order; filled XFA form datasets or static template text when
the PDF only shows a "Please wait" page), .docx (paragraphs; table cells
tab-separated), .xlsx (every sheet, tab-separated, under '## sheet <name>'), .html/.htm
(tags stripped, block elements become line breaks), .xml (JATS full text from PMC /
Europe PMC / eLife: tags stripped, paragraphs/titles/table rows on their own lines, cells
tab-separated), and -- when LibreOffice (`soffice`) is installed -- the legacy binary
.xls/.doc/.ppt and .rtf documents (converted to xlsx/docx/pdf first, then read as those). A PDF with no text
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
import xml.etree.ElementTree as ET
import zlib
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from scrape_primers import HTML_EXT, OFFICE_EXT  # noqa: E402
from scrape_primers import read as _read  # noqa: E402

XML_EXT = {".xml"}
LEGACY_EXT = {".xls", ".doc", ".ppt", ".rtf"}
SOFFICE = shutil.which("soffice") or shutil.which("libreoffice")
GHOSTSCRIPT = shutil.which("gs")
CONVERT = {".pdf", ".eps"} | OFFICE_EXT | HTML_EXT | XML_EXT | (LEGACY_EXT if SOFFICE else set())

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
    """Office/RTF via LibreOffice -> xlsx/docx/pdf in a scratch dir, then read that."""
    if not SOFFICE:
        print(f"skip {path}: legacy {path.suffix} needs LibreOffice (soffice)", file=sys.stderr)
        return None
    target = {".xls": "xlsx", ".doc": "docx", ".ppt": "pdf", ".rtf": "docx"}[path.suffix.lower()]
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


def _eps(path: Path) -> str | None:
    """Read an EPS figure through an optional local Ghostscript conversion."""
    if not GHOSTSCRIPT:
        print(f"skip {path}: EPS text extraction needs Ghostscript (gs)", file=sys.stderr)
        return None
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "figure.pdf"
        try:
            done = subprocess.run([GHOSTSCRIPT, "-dSAFER", "-dBATCH", "-dNOPAUSE",
                                   "-dEPSCrop", "-sDEVICE=pdfwrite",
                                   f"-sOutputFile={out}", str(path.resolve())],
                                  capture_output=True, timeout=90, check=False)
        except (OSError, subprocess.TimeoutExpired):
            print(f"skip {path}: EPS conversion did not complete", file=sys.stderr)
            return None
        if done.returncode or not out.is_file() or out.stat().st_size == 0:
            print(f"skip {path}: Ghostscript could not convert EPS", file=sys.stderr)
            return None
        return read(out)


def xfa_text(raw: bytes) -> str | None:
    """Read filled Adobe LiveCycle form data from Flate-compressed PDF streams.

    Ordinary PDF text extraction often sees only the static "Please wait" page. XFA
    may hold both a blank form template and a later filled dataset, so select the
    dataset with the most nonempty field text. If there is no filled dataset, retain
    static text values from the form template, explicitly labelled as such. This
    reads data only; form scripts and presentation markup are not executed.
    """
    datasets: list[tuple[int, str]] = []
    templates: list[str] = []
    for match in re.finditer(rb"\b\d+\s+\d+\s+obj\s*<<(.*?)>>\s*stream\r?\n", raw, re.S):
        end = raw.find(b"endstream", match.end())
        if end < 0:
            continue
        body = raw[match.end():end].rstrip(b"\r\n")
        try:
            decoded = zlib.decompress(body) if b"/FlateDecode" in match.group(1) else body
        except zlib.error:
            continue
        is_dataset = b"datasets" in decoded[:256] and b"xfa-data" in decoded[:256]
        is_template = b"xfa-template" in decoded[:256]
        if not (is_dataset or is_template):
            continue
        try:
            tree = ET.fromstring(decoded)
        except ET.ParseError:
            continue
        if is_template:
            values = []
            for node in tree.iter():
                if node.tag.rsplit("}", 1)[-1] != "value":
                    continue
                for value in node:
                    tag = value.tag.rsplit("}", 1)[-1]
                    if tag not in {"text", "exData"}:
                        continue
                    if tag == "exData" and value.get("contentType", "").lower() not in {"text/html", "text/plain"}:
                        continue
                    text = "".join(value.itertext()).strip()
                    if text:
                        values.append(text)
            if values:
                templates.append("\n".join(values) + "\n")
            continue
        fields: list[str] = []

        def walk(node: ET.Element, path: tuple[str, ...]) -> None:
            tag = node.tag.rsplit("}", 1)[-1]
            if tag not in {"datasets", "data", "form"}:
                path += (tag,)
            if len(node) == 0 and node.text and node.text.strip():
                value = node.text.strip().replace("\r", "\n")
                fields.append("/".join(path) + ": " + value)
            for child in node:
                walk(child, path)

        walk(tree, ())
        if fields:
            datasets.append((sum(len(f) for f in fields), "\n".join(fields) + "\n"))
    if datasets:
        return max(datasets)[1]
    if templates:
        return ("## XFA template text (form layout and scripts not evaluated)\n"
                + max(templates, key=len))
    return None


def read(path: Path) -> str | None:
    """Text of a document (None if unreadable). Extends scrape_primers.read with XML and
    the legacy Office formats."""
    ext = path.suffix.lower()
    if ext in XML_EXT:
        return jats_text(path.read_text(encoding="utf-8", errors="replace"))
    if ext in LEGACY_EXT:
        return _legacy(path)
    if ext == ".eps":
        return _eps(path)
    text = _read(path)
    if ext == ".pdf" and (text is None or len(text.strip()) < 100
                          or text.lstrip().startswith("Please wait...")):
        xfa = xfa_text(path.read_bytes())
        if xfa is not None:
            return xfa
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
