"""Plumbing shared by get_sources.py and every fetcher: HTTP (curl), download validation,
and the per-protocol source Folder with its MANIFEST.tsv.

Nothing here knows about a particular publisher. The one rule it enforces is that a file
is only kept when it is what its name says: a .pdf starts with %PDF and ends with %%EOF,
an .xlsx/.docx/.zip is a complete zip archive, an .html page is not a captcha / Cloudflare
/ "Preparing to download" interstitial and has real text, and so on (`validate`). Anything
else is recorded in the manifest as `(manual) NAME` with the URL and the reason, so a
person (or a later run) can fetch it by hand.
"""
from __future__ import annotations

import gzip
import hashlib
import html
from html.parser import HTMLParser
import json
import re
import shutil
import subprocess
import sys
import time
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "lib"))
from storage import data_root

UA = "Mozilla/5.0 (X11; Linux x86_64) chem-sources/1.0 (research; mailto:he.johan@gmail.com)"
MAX_BYTES = 120_000_000     # per file; bigger is a data dump (counts, peaks), not a methods source
CONVERT_MAX = 15_000_000    # spreadsheets/zips above this are not turned into .txt automatically
RETRY_STATUS = {0, 429, 500, 502, 503, 504}
LEGACY_CAPS = {50_000_000}  # the old curl --max-filesize: a file of exactly this size was cut off

# Interstitials that come back with HTTP 200. Matched against the first 64 kB.
GATE = re.compile(
    rb"<title>\s*(?:Just a moment|Checking your browser|Preparing to download|Access Denied"
    rb"|Attention Required|Are you a robot|Security check|Please wait)"
    rb"|recaptcha/challengepage|RecaptchaChallengePageUi|cf-browser-verification"
    rb"|/cdn-cgi/challenge-platform|cf_chl_opt|Enable JavaScript and cookies to continue"
    rb"|Checking your browser before accessing", re.I)

PDF_EXT = {".pdf"}
ZIP_EXT = {".zip", ".xlsx", ".xlsm", ".docx", ".pptx"}
OLE_EXT = {".xls", ".doc", ".ppt"}
HTML_EXT = {".html", ".htm"}
XML_EXT = {".xml"}
JSON_EXT = {".json"}
GZ_EXT = {".gz", ".tgz"}
TEXT_EXT = {".txt", ".csv", ".tsv", ".fa", ".fasta", ".fastq", ".md", ".tab", ".bed", ".gb"}
MEDIA_EXT = {".png", ".jpg", ".jpeg", ".gif", ".tif", ".tiff", ".avi", ".mp4", ".mov",
             ".dwg", ".dxf", ".svg", ".eps"}


def say(*a) -> None:
    print(*a, flush=True)


def data_dir() -> Path:
    return data_root() / "sources"


def cache_dir() -> Path:
    d = data_dir() / ".cache"
    d.mkdir(parents=True, exist_ok=True)
    return d


# ------------------------------------------------------------------------- http

def get(url: str, pause: float = 0.4, headers: dict | None = None,
        max_time: int = 120) -> tuple[int, bytes]:
    """-> (HTTP status, body) for small API calls. curl: the system Python has no CA bundle."""
    time.sleep(pause)
    cmd = ["curl", "-sSL", "--connect-timeout", "30", "--max-time", str(max_time),
           "--max-filesize", str(MAX_BYTES), "-A", UA, "-w", "\n%{http_code}"]
    for k, v in (headers or {}).items():
        cmd += ["-H", f"{k}: {v}"]
    done = subprocess.run(cmd + [url], capture_output=True, check=False)
    body, _, code = done.stdout.rpartition(b"\n")
    return (int(code) if code.isdigit() else 0), body


def get_json(url: str, headers: dict | None = None, retries: int = 1) -> dict:
    for attempt in range(retries + 1):
        code, body = get(url, headers=headers)
        if code == 200:
            try:
                return json.loads(body)
            except ValueError:
                return {}
        if code not in RETRY_STATUS or attempt == retries:
            return {}
        time.sleep(3 * (attempt + 1))
    return {}


def download(url: str, dest: Path, max_bytes: int = MAX_BYTES,
             max_time: int = 900) -> tuple[int, str]:
    """curl `url` into `dest`. -> (HTTP status, problem or "")."""
    time.sleep(0.4)
    done = subprocess.run(["curl", "-sSL", "--connect-timeout", "30", "--max-time",
                           str(max_time), "--max-filesize", str(max_bytes), "-A", UA,
                           "-o", str(dest), "-w", "%{http_code}", url],
                          capture_output=True, check=False)
    out = done.stdout.decode(errors="replace").strip()
    code = int(out[-3:]) if out[-3:].isdigit() else 0
    problem = {0: "", 63: f"larger than {max_bytes // 1_000_000} MB cap (not fetched)",
               28: "timed out (download incomplete)", 18: "transfer closed early (truncated)",
               }.get(done.returncode, f"curl error {done.returncode}")
    return code, problem


# ------------------------------------------------------------------- validation

def visible_text(raw: bytes) -> str:
    t = raw.decode("utf-8", errors="replace")
    t = re.sub(r"(?is)<(script|style|noscript|svg)\b.*?</\1>", " ", t)
    t = re.sub(r"<[^>]+>", " ", t)
    return re.sub(r"\s+", " ", html.unescape(t)).strip()


def madcap_topic_text(raw: bytes) -> str:
    """Recognize a complete short documentation topic, excluding its navigation.

    An HTML5 help-system shell is not the chapter it loads. Only a topic document
    with a closed main-content container and substantive text qualifies here.
    """
    class Topic(HTMLParser):
        def __init__(self):
            super().__init__()
            self.topic = False
            self.depth = 0
            self.started = False
            self.hidden = 0
            self.parts = []

        def handle_starttag(self, tag, attrs):
            attrs = dict(attrs)
            if tag == "html":
                self.topic = "Topic" in attrs.get("data-mc-runtime-file-type", "").split(";")
            if tag == "div":
                if self.depth:
                    self.depth += 1
                elif not self.started and attrs.get("id") == "mc-main-content":
                    self.started = True
                    self.depth = 1
            if self.depth and tag in {"script", "style", "noscript", "svg"}:
                self.hidden += 1

        def handle_endtag(self, tag):
            if self.depth and tag in {"script", "style", "noscript", "svg"}:
                self.hidden = max(0, self.hidden - 1)
            if tag == "div" and self.depth:
                self.depth -= 1

        def handle_data(self, data):
            if self.depth and not self.hidden:
                self.parts.append(data)

    parser = Topic()
    parser.feed(raw.decode("utf-8", errors="replace"))
    text = re.sub(r"\s+", " ", " ".join(parser.parts)).strip()
    return text if parser.topic and parser.started and parser.depth == 0 and len(text) >= 100 else ""


def looks_html(head: bytes) -> bool:
    h = head.lstrip()[:1000].lower()
    return h.startswith((b"<!doctype html", b"<html")) or (b"<html" in h and b"<head" in h)


def gate_reason(head: bytes) -> str:
    m = GATE.search(head)
    if not m:
        return ""
    s = m.group(0).decode(errors="replace")
    if "recaptcha" in s.lower():
        return "reCAPTCHA challenge page"
    if "Preparing to download" in s:
        return "PMC 'Preparing to download' gate"
    return f"challenge/gate page ({re.sub(r'<[^>]+>', '', s).strip()[:40]})"


def validate(path: Path, name: str | None = None, min_text: int = 1500,
             need: bytes | None = None) -> str:
    """'' if the file at `path` is a genuine instance of what `name`'s extension promises,
    else the reason it is not. `need`: a byte regex that must occur (e.g. rb'<body' for a
    full-text JATS XML). `min_text`: least visible text an HTML page must carry."""
    name = name or path.name
    ext = Path(name.lower().removesuffix(".part")).suffix
    try:
        size = path.stat().st_size
    except OSError:
        return "missing"
    if size == 0:
        return "empty file"
    if size in LEGACY_CAPS:
        return f"exactly {size:,} bytes: cut off at an older get_sources size cap (truncated)"
    with path.open("rb") as fh:
        head = fh.read(65536)
        fh.seek(max(0, size - 32768))
        tail = fh.read()
    gate = gate_reason(head)
    if ext in PDF_EXT:
        if head.lstrip()[:5] != b"%PDF-" and b"%PDF-" not in head[:1024]:
            return gate or ("HTML page, not a PDF" if looks_html(head) else "not a PDF (no %PDF header)")
        if b"%%EOF" not in tail:
            return "truncated PDF (no %%EOF trailer)"
        return ""
    if ext in ZIP_EXT:
        if head[:2] != b"PK":
            return gate or ("HTML page, not a zip/Office file" if looks_html(head)
                            else "not a zip/Office file")
        try:
            with zipfile.ZipFile(path) as z:
                bad = z.testzip()
                if bad:
                    return f"corrupt zip member {bad}"
                if ext == ".zip" and not [i for i in z.infolist() if i.file_size]:
                    return "zip holds no files"
        except (zipfile.BadZipFile, OSError, EOFError, zipfile.LargeZipFile) as e:
            return f"truncated or corrupt zip ({e})"
        return ""
    if ext in OLE_EXT:
        if head[:8] == bytes.fromhex("D0CF11E0A1B11AE1") or head[:2] == b"PK":
            return ""
        return gate or ("HTML page, not a legacy Office file" if looks_html(head)
                        else "not a legacy Office (OLE2) file")
    if ext in GZ_EXT:
        if head[:2] != b"\x1f\x8b":
            return gate or "not gzip data"
        try:
            with gzip.open(path, "rb") as g:
                while g.read(1 << 20):
                    pass
        except (EOFError, OSError, gzip.BadGzipFile) as e:
            return f"truncated or corrupt gzip ({e})"
        return ""
    if gate:
        return gate
    if ext in HTML_EXT:
        raw = head if size <= 65536 else path.read_bytes()
        n = len(visible_text(raw))
        if n < min_text and not madcap_topic_text(raw):
            return f"no article content ({n} characters of text)"
        if need and not re.search(need, path.read_bytes()):
            return f"expected content missing ({need.decode(errors='replace')})"
        return ""
    if ext in XML_EXT:
        s = head.lstrip()
        if not s.startswith(b"<"):
            return "not XML"
        if b"errorBean" in head[:2000] or b"<Error>" in head[:400]:
            return "API error document, not the content"
        if looks_html(head):
            return "HTML page, not XML"
        if need and not re.search(need, path.read_bytes()):
            return f"expected content missing ({need.decode(errors='replace')}): front matter only"
        return ""
    if ext in JSON_EXT:
        try:
            json.loads(path.read_bytes())
        except ValueError:
            return "HTML page, not JSON" if looks_html(head) else "not valid JSON"
        return ""
    if ext in TEXT_EXT:
        if b"\0" in head:
            return "binary data, not text"
        if looks_html(head):
            return "HTML page, not a text file"
        return ""
    if ext in MEDIA_EXT and looks_html(head):
        return "HTML page, not a media file"
    return ""


# ------------------------------------------------------------------------ folder

def safe_name(name: str) -> str:
    return re.sub(r"[^\w.+-]", "_", name)[:140]


def md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def is_twin(path: Path) -> bool:
    """FILE.ext.txt written by doctext beside FILE.ext (or beside its demoted copy)."""
    if not path.name.endswith(".txt"):
        return False
    stem = path.with_name(path.name[:-4])
    return stem.exists() or stem.suffix.lower() in (PDF_EXT | ZIP_EXT | OLE_EXT | HTML_EXT
                                                     | XML_EXT | {".txt"})


class Folder:
    """One protocol's source directory and its manifest (file, url, bytes, md5, what)."""

    HEADER = "file\turl\tbytes\tmd5\twhat\n"

    def __init__(self, path: Path, convert: bool = True, expand: bool = True):
        self.path = path
        self.convert = convert
        self.expand = expand        # unpack the documents inside a fetched .zip
        path.mkdir(parents=True, exist_ok=True)
        self.manifest = path / "MANIFEST.tsv"
        self.rows: dict[str, list[str]] = {}
        self.gated_hosts: set[str] = set()
        if self.manifest.exists():
            for line in self.manifest.read_text(errors="replace").splitlines()[1:]:
                f = line.split("\t")
                if f and f[0]:
                    self.rows[f[0]] = (f + [""] * 5)[:5]
        self.by_md5 = {r[3]: r[0] for r in self.rows.values()
                       if r[3] and not r[0].startswith("(")}

    # --- bookkeeping
    def have(self, name: str) -> bool:
        return (self.path / name).exists() and name in self.rows

    def manual(self, name: str, url: str, what: str, reason: str) -> None:
        if name in self.rows and self.rows[name][1] == url and (self.path / name).exists():
            return
        key, alternate = self._manual_keys(name, url)
        if ((key in self.rows and self.rows[key][1] != url) or alternate in self.rows
                or (name in self.rows and self.rows[name][1] != url)):
            key = alternate
        self.rows[key] = [key, url, "", "", f"{what}; {reason}"]
        say(f"    manual  {name}  ({reason})  <- {url}")
        self.write()

    @staticmethod
    def _manual_keys(name: str, url: str) -> tuple[str, str]:
        path = Path(name)
        tagged = path.stem + "__" + hashlib.sha256(url.encode()).hexdigest()[:12] + path.suffix
        return f"(manual) {name}", f"(manual) {tagged}"

    def _clear_manual(self, name: str, url: str) -> None:
        for key in self._manual_keys(name, url):
            if key in self.rows and self.rows[key][1] == url:
                self.rows.pop(key)

    def note(self, key: str, url: str, what: str) -> None:
        self.rows[f"(note) {key}"] = [f"(note) {key}", url, "", "", what]

    def _record(self, name: str, url: str, what: str) -> bool:
        """Row for a validated file now at self.path/name. False if it duplicates one."""
        p = self.path / name
        digest = md5(p)
        twin = self.by_md5.get(digest)
        if twin and twin != name and (self.path / twin).exists():
            p.unlink()
            self.rows[f"(duplicate) {name}"] = [f"(duplicate) {name}", url, "", digest,
                                                f"{what}; identical to {twin}"]
            self._clear_manual(name, url)
            say(f"    dup     {name}  = {twin}")
            self.write()
            return False
        self.by_md5[digest] = name
        self._clear_manual(name, url)
        self.rows[name] = [name, url, str(p.stat().st_size), digest, what]
        self.write()
        return True

    # --- fetching
    def _source_collision(self, name: str, url: str, what: str) -> bool:
        """A valid file at the same name does not establish source identity."""
        row = self.rows.get(name) or self.rows.get(f"(duplicate) {name}")
        if row is None or not row[1] or row[1] == url:
            return False
        key = "source-collision-" + name + "-" + hashlib.sha256(url.encode()).hexdigest()[:12]
        self.manual(key, url, what,
                    f"cache filename {name} is already assigned to {row[1]}; "
                    "use a distinct source filename; existing file was not reused or overwritten")
        return True

    def save(self, name: str, url: str, what: str, *, min_text: int = 1500,
             need: bytes | None = None, max_bytes: int = MAX_BYTES,
             retries: int = 2) -> Path | None:
        """Download `url` as `name` unless already here; validate; convert to .txt.
        -> the saved path, or None (then a `(manual)` row says why)."""
        name = safe_name(name)
        if self._source_collision(name, url, what):
            return None
        dest = self.path / name
        if self.have(name) and not validate(dest, name, min_text, need):
            return dest
        if dest.exists() and not validate(dest, name, min_text, need):  # on disk, no row
            if self._record(name, url, what + " (already on disk)"):
                self.to_text(dest)
                return dest
            return None
        if f"(duplicate) {name}" in self.rows:
            return None
        host = re.sub(r"^\w+://([^/]+).*", r"\1", url)
        part = self.path / (name + ".part")
        reason, code = "", 0
        for attempt in range(retries + 1):
            code, problem = download(url, part, max_bytes)
            if problem:
                reason = problem
            elif code != 200:
                reason = f"HTTP {code}"
            else:
                reason = validate(part, name, min_text, need)
            if not reason:
                break
            gated = "gate" in reason or "challenge" in reason or "reCAPTCHA" in reason
            if attempt == retries or not (gated or code in RETRY_STATUS) or "cap" in reason:
                break
            time.sleep(5 * (attempt + 1) ** 2)
        if reason:
            part.unlink(missing_ok=True)
            if "gate" in reason or "challenge" in reason or "reCAPTCHA" in reason:
                self.gated_hosts.add(host)
            self.manual(name, url, what, reason)
            return None
        part.replace(dest)
        if not self._record(name, url, what):
            return None
        say(f"    saved   {name}  ({dest.stat().st_size:,} bytes)")
        self.to_text(dest)
        if dest.suffix.lower() == ".zip" and self.expand:
            self.expand_zip(name, dest.stem + "_", url, what)
        return dest

    def keep(self, name: str, data: bytes, url: str, what: str) -> Path | None:
        """Store bytes we already hold (from a zip, a cache) -- validated like a download."""
        name = safe_name(name)
        if self._source_collision(name, url, what):
            return None
        dest = self.path / name
        part = self.path / (name + ".part")
        part.write_bytes(data)
        reason = validate(part, name)
        if reason:
            part.unlink()
            self.manual(name, url, what, reason)
            return None
        part.replace(dest)
        if not self._record(name, url, what):
            return None
        say(f"    saved   {name}  ({len(data):,} bytes)")
        self.to_text(dest)
        return dest

    def to_text(self, src: Path) -> None:
        """The plain-text twin, right away, so a reader can start before the run ends."""
        if not self.convert:
            return
        try:
            import doctext
        except ImportError:
            return
        if src.suffix.lower() not in doctext.CONVERT:
            return
        size = src.stat().st_size
        if size > CONVERT_MAX and src.suffix.lower() not in {".pdf", ".html", ".htm", ".xml"}:
            row = self.rows.get(src.name)
            if row and "not converted" not in row[4]:
                row[4] += (f"; not converted to text ({size // 1_000_000} MB: run "
                           f"tools/doctext.py on it if needed)")
                self.write()
            say(f"    (no .txt: {size // 1_000_000} MB {src.suffix})")
            return
        try:
            out = doctext.convert(src)
        except Exception as e:  # noqa: BLE001 -- a bad document must not stop the run
            say(f"    (no .txt: {e})")
            return
        if out:
            say(f"    text    {out.name}  ({out.stat().st_size:,} bytes)")

    def expand_zip(self, name: str, prefix: str, url: str, what: str) -> None:
        """Pull readable documents out of a saved zip into top-level, tracked files."""
        keep_ext = (PDF_EXT | ZIP_EXT | OLE_EXT | HTML_EXT | XML_EXT | TEXT_EXT) - {".zip"}
        try:
            with zipfile.ZipFile(self.path / name) as z:
                for info in z.infolist():
                    m = info.filename
                    base = Path(m).name
                    if (m.endswith("/") or not info.file_size or base.startswith((".", "__"))
                            or "__MACOSX" in m or Path(base).suffix.lower() not in keep_ext
                            or info.file_size > 40_000_000):
                        continue
                    target = safe_name(prefix + base)
                    if (self.path / target).exists() or f"(duplicate) {target}" in self.rows:
                        continue
                    self.keep(target, z.read(m), f"{url}#{m}", f"{what} (from the zip)")
        except (zipfile.BadZipFile, OSError) as e:
            say(f"    (cannot open {name}: {e})")

    def write(self) -> None:
        tmp = self.manifest.with_suffix(".tsv.part")
        with tmp.open("w") as fh:
            fh.write(self.HEADER)
            for r in sorted(self.rows.values()):
                fh.write("\t".join(c.replace("\t", " ").replace("\n", " ") for c in r) + "\n")
        tmp.replace(self.manifest)


def have_tool(name: str) -> bool:
    return shutil.which(name) is not None
