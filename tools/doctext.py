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
.xls/.doc/.ppt and .rtf documents (converted to xlsx/docx/pdf first, then read as those;
a .doc falls back to catdoc). An OLE2 legacy file saved under an OOXML name (.docx/.xlsx)
is recognised by its bytes. Huge sheets are streamed. A PDF with no text layer (a scanned
or "Print to PDF" image) is OCR'd when tesseract and pdftoppm are installed: the twin
then starts with an "## OCR TEXT" header; a single image-only page inside a text PDF
(a pasted-in table) is OCR'd and appended after the text layer under its own header.
OCR text is a reading aid: verify it against the rendered original. A .txt is rewritten
only when its source is newer. tools/get_sources.py runs this on everything it downloads.
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
CATDOC = shutil.which("catdoc")
GHOSTSCRIPT = shutil.which("gs")
TESSERACT = shutil.which("tesseract")
PDFTOPPM = shutil.which("pdftoppm")
PDFIMAGES = shutil.which("pdfimages")
PDFINFO = shutil.which("pdfinfo")
CONVERT = ({".pdf", ".eps"} | OFFICE_EXT | HTML_EXT | XML_EXT | (LEGACY_EXT if SOFFICE else set())
           | ({".doc"} if CATDOC else set()))
OLE_MAGIC = bytes.fromhex("D0CF11E0A1B11AE1")
OLE_FOR_OOXML = {".docx": ".doc", ".xlsx": ".xls", ".xlsm": ".xls", ".pptx": ".ppt"}

OCR_EMPTY = 20             # fewer non-blank characters than this: the PDF has no text layer
OCR_MAX_PAGES = 40          # a longer image-only PDF is OCR'd up to here, with a note
OCR_DPI = 300
OCR_MAX_PIXELS = 4200       # longest rendered side: a poster-sized page is rendered coarser
OCR_HEADER = ("## OCR TEXT (tesseract, rendered up to {dpi} dpi): machine-read from the rendered pages, not a "
              "text layer -- verify every base, number and symbol against the rendered original")
OCR_PAGE_HEADER = ("## OCR TEXT of image-only page {page} (tesseract, up to {dpi} dpi): machine-read, "
                   "appended after the text layer -- verify against the rendered page")

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


def _legacy(path: Path, suffix: str = "") -> str | None:
    """Office/RTF via LibreOffice -> xlsx/docx/pdf in a scratch dir, then read that. A
    legacy Word .doc falls back to catdoc. `suffix` overrides a misleading file suffix
    (an OLE2 .doc saved under a .docx name)."""
    renamed = bool(suffix) and suffix != path.suffix.lower()
    suffix = suffix or path.suffix.lower()
    target = {".xls": "xlsx", ".doc": "docx", ".ppt": "pdf", ".rtf": "docx"}[suffix]
    if SOFFICE:
        with tempfile.TemporaryDirectory() as tmp:
            src = path.resolve()
            if renamed:                         # give LibreOffice the true suffix
                src = Path(tmp) / "in" / (path.stem + suffix)
                src.parent.mkdir()
                shutil.copyfile(path, src)
            profile = Path(tmp) / "profile"      # private profile: parallel runs do not collide
            try:
                done = subprocess.run([SOFFICE, f"-env:UserInstallation=file://{profile}", "--headless",
                                       "--convert-to", target, "--outdir", tmp, str(src)],
                                      capture_output=True, timeout=300, check=False)
            except (OSError, subprocess.TimeoutExpired):
                done = None
            out = Path(tmp) / (path.stem + "." + target)
            if done is not None and not done.returncode and out.exists():
                return _read(out)
    if suffix == ".doc" and CATDOC:
        try:
            done = subprocess.run([CATDOC, "-w", "-d", "utf-8", str(path)], capture_output=True,
                                  timeout=120, check=False)
        except (OSError, subprocess.TimeoutExpired):
            done = None
        if done is not None and not done.returncode and done.stdout.strip():
            return done.stdout.decode("utf-8", errors="replace")
    need = "LibreOffice (soffice)" + (" or catdoc" if suffix == ".doc" else "")
    print(f"skip {path}: legacy {suffix} needs {need}, or it could not be converted",
          file=sys.stderr)
    return None


def _ole_suffix(path: Path) -> str:
    """The legacy suffix of an OLE2 Office file saved under an OOXML name, else ''."""
    legacy = OLE_FOR_OOXML.get(path.suffix.lower())
    if not legacy:
        return ""
    with path.open("rb") as fh:
        return legacy if fh.read(8) == OLE_MAGIC else ""


_ROW = re.compile(r"<row\b[^>]*?(?:/>|>(.*?)</row>)", re.S)
_CELL = re.compile(r"<c\b([^>]*?)(?:/>|>(.*?)</c>)", re.S)


def xlsx_text(path: Path) -> str:
    """Every sheet as tab-separated rows under '## sheet <name>' (the format of
    scrape_primers.xlsx_text), streaming each sheet's XML in chunks: a sheet of a million
    styled but empty rows (self-closing <row/>) is skipped row by row instead of being
    read and scanned as one string."""
    import zipfile
    from scrape_primers import _xml_text
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

        def row_line(row: str) -> None:
            cells = []
            for attrs, body in _CELL.findall(row):
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

        for i, n in enumerate(sheets):
            out.append(f"## sheet {html.unescape(titles[i]) if i < len(titles) else i + 1}")
            buf = ""
            with z.open(n) as fh:
                while True:
                    chunk = fh.read(1 << 20)
                    buf += chunk.decode("utf-8", "replace") if chunk else ""
                    end = 0
                    for m in _ROW.finditer(buf):
                        if m.group(1):
                            row_line(m.group(1))
                        end = m.end()
                    # keep only an unfinished row (from its last '<row' start) for the next chunk
                    rest = buf[end:]
                    cut = rest.rfind("<row")
                    buf = rest[cut:] if cut >= 0 else rest[-16:]
                    if not chunk:
                        break
        return "\n".join(out) + "\n"


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
        return read(out, ocr=False)      # a figure: its own text only, not OCR of its art


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


def _run(args: list[str], timeout: int) -> subprocess.CompletedProcess | None:
    try:
        return subprocess.run(args, capture_output=True, timeout=timeout, check=False)
    except (OSError, subprocess.TimeoutExpired):
        return None


def ocr_page(path: Path, page: int, longest_pt: float = 792) -> str | None:
    """OCR one rendered PDF page (pdftoppm + tesseract); None when unavailable or failed.
    Rendered at OCR_DPI, or coarser so the longest side stays within OCR_MAX_PIXELS (a
    2560x3300 pt page holding a 72 ppi scan has nothing more to show at 300 dpi)."""
    if not (TESSERACT and PDFTOPPM):
        return None
    dpi = max(50, min(OCR_DPI, int(OCR_MAX_PIXELS * 72 / max(longest_pt, 1))))
    with tempfile.TemporaryDirectory() as tmp:
        done = _run([PDFTOPPM, "-f", str(page), "-l", str(page), "-r", str(dpi), "-gray",
                     "-png", str(path.resolve()), str(Path(tmp) / "page")], 180)
        images = sorted(Path(tmp).glob("page*.png"))
        if done is None or done.returncode or not images:
            return None
        done = _run([TESSERACT, str(images[0]), "-"], 300)
        if done is None or done.returncode:
            return None
        return done.stdout.decode("utf-8", errors="replace")


def _page_sizes(path: Path, pages: int) -> dict[int, tuple[float, float]]:
    """{page: (width, height) in points} from pdfinfo."""
    done = _run([PDFINFO, "-f", "1", "-l", str(pages), str(path)], 60) if PDFINFO else None
    out = done.stdout.decode(errors="replace") if done is not None else ""
    return {int(p): (float(w), float(h)) for p, w, h in
            re.findall(r"Page\s+(\d+)\s+size:\s+([\d.]+)\s+x\s+([\d.]+)", out)}


def _longest(sizes: dict, page: int) -> float:
    return max(sizes.get(page, (612.0, 792.0)))


def raster_pages(path: Path, pages: list[int]) -> list[int]:
    """Those of `pages` whose content is one raster image covering at least half the page
    (a scanned or pasted-in page), from `pdfimages -list`."""
    if not (PDFIMAGES and pages):
        return []
    done = _run([PDFIMAGES, "-list", "-f", str(min(pages)), "-l", str(max(pages)), str(path)], 120)
    if done is None or done.returncode:
        return []
    sizes = _page_sizes(path, max(pages))
    found = set()
    for line in done.stdout.decode(errors="replace").splitlines()[2:]:
        f = line.split()
        try:
            page, width, height, xppi, yppi = int(f[0]), int(f[3]), int(f[4]), int(f[12]), int(f[13])
        except (IndexError, ValueError):
            continue
        w, h = sizes.get(page, (0.0, 0.0))
        area = w * h / 72 / 72
        if page in pages and area and xppi and yppi and width / xppi * height / yppi >= 0.5 * area:
            found.add(page)
    return sorted(found)


def with_ocr(path: Path, text: str) -> str:
    """Add OCR text where a PDF has no text layer. A PDF without text (fewer than OCR_EMPTY
    non-blank characters, e.g. only form feeds) is OCR'd page by
    page (up to OCR_MAX_PAGES) and the twin starts with an OCR header. In a PDF with a text
    layer, a page that is one raster image with almost no text (a table pasted in as a
    picture) is OCR'd and appended after the text layer, under its own header, so line
    numbers of the text layer stay stable."""
    if not (TESSERACT and PDFTOPPM):
        if len(re.sub(r"\s", "", text)) < OCR_EMPTY:
            print(f"warning {path}: no text layer (image-only PDF?) and no OCR tools "
                  f"(tesseract, pdftoppm) -- read the rendered pages", file=sys.stderr)
        return text
    pages = text.split("\f")
    if pages and not pages[-1].strip():
        pages = pages[:-1]
    if len(re.sub(r"\s", "", text)) < OCR_EMPTY:
        n = len(pages) or 1
        sizes = _page_sizes(path, min(n, OCR_MAX_PAGES))
        done = [ocr_page(path, p, _longest(sizes, p)) for p in range(1, min(n, OCR_MAX_PAGES) + 1)]
        if not any(t and t.strip() for t in done):
            print(f"warning {path}: no text layer and OCR found no text", file=sys.stderr)
            return text
        head = OCR_HEADER.format(dpi=OCR_DPI)
        if n > OCR_MAX_PAGES:
            head += f"\n## only pages 1-{OCR_MAX_PAGES} of {n} OCR'd (doctext OCR_MAX_PAGES)"
        return head + "\n\n" + "\f".join(t or "" for t in done)
    sparse = [i + 1 for i, t in enumerate(pages) if len(re.sub(r"\s", "", t)) < 400]
    extra = []
    raster = raster_pages(path, sparse)[:OCR_MAX_PAGES]
    sizes = _page_sizes(path, max(raster)) if raster else {}
    for page in raster:
        t = ocr_page(path, page, _longest(sizes, page))
        if t and t.strip():
            extra.append(OCR_PAGE_HEADER.format(page=page, dpi=OCR_DPI) + "\n\n" + t.strip() + "\n")
    return text if not extra else text.rstrip("\n") + "\n\n" + "\n".join(extra)


def read(path: Path, ocr: bool = True) -> str | None:
    """Text of a document (None if unreadable). Extends scrape_primers.read with XML, the
    legacy Office formats (also under a misleading OOXML suffix), streamed XLSX and OCR of
    PDF pages without a text layer."""
    ext = path.suffix.lower()
    if ext in XML_EXT:
        return jats_text(path.read_text(encoding="utf-8", errors="replace"))
    if ext in LEGACY_EXT:
        return _legacy(path)
    if ext == ".eps":
        return _eps(path)
    legacy = _ole_suffix(path)
    if legacy:
        return _legacy(path, legacy)
    if ext in {".xlsx", ".xlsm"}:
        try:
            return xlsx_text(path)
        except Exception as e:                          # noqa: BLE001
            print(f"skip {path}: cannot read ({e})", file=sys.stderr)
            return None
    text = _read(path)
    if ext == ".pdf" and (text is None or len(text.strip()) < 100
                          or text.lstrip().startswith("Please wait...")):
        xfa = xfa_text(path.read_bytes())
        if xfa is not None:
            return xfa
    if text is not None and ext == ".pdf" and ocr:
        text = with_ocr(path, text)
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
