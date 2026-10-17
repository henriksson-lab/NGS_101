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
import fcntl
import subprocess
import sys
import time
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "tools"))
sys.path.insert(0, str(ROOT / "lib"))
from storage import data_root

UA = "Mozilla/5.0 (X11; Linux x86_64) chem-sources/1.0 (research)"
MAX_BYTES = 120_000_000     # per file; bigger is a data dump (counts, peaks), not a methods source
CONVERT_MAX = 15_000_000    # spreadsheets/zips above this are not turned into .txt automatically
RETRY_STATUS = {0, 429, 500, 502, 503, 504}
LEGACY_CAPS = {50_000_000}  # the old curl --max-filesize: a file of exactly this size was cut off

# Interstitials that come back with HTTP 200. Matched against the first 64 kB.
GATE = re.compile(
    rb"<title>\s*(?:Just a moment|Checking your browser|Preparing to download|Access Denied|Client Challenge"
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
OLE_MAGIC = bytes.fromhex("D0CF11E0A1B11AE1")
# An OOXML name whose bytes are the legacy OLE2 container: stored under this suffix.
OLE_FOR_OOXML = {".docx": ".doc", ".xlsx": ".xls", ".pptx": ".ppt", ".xlsm": ".xls"}
# Error bodies served with HTTP 200 that are no document of any kind.
ERROR_BODY = re.compile(rb"PDF link has expired or is invalid", re.I)
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

STALL_SECONDS = 60          # a download that receives nothing for this long is abandoned
CURL_HTTP2_STREAM = 92     # curl: HTTP/2 stream error (some vendor CDNs); retried once as HTTP/1.1


# Polite access: at most one request per host per interval, shared by every fetcher process
# through lock files under the cache (parallel runs used to burst a publisher into serving a
# challenge). A host that serves a challenge or 429 is left alone for a cool-down period instead
# of being retried; a challenge page is an access decision, never something to work around.
HOST_INTERVAL = 1.0
SLOW_HOSTS = ("springer.com", "nature.com", "springernature.com", "wiley.com", "sciencedirect.com",
              "els-cdn.com", "cell.com", "science.org", "pnas.org", "oup.com", "cshlp.org",
              "biorxiv.org", "medrxiv.org", "researchsquare.com", "tandfonline.com", "asm.org")
SLOW_INTERVAL = 8.0
COOLDOWN_SECONDS = 6 * 3600


def host_of(url: str) -> str:
    return re.sub(r"^\w+://([^/:]+).*", r"\1", url).lower()


def host_interval(host: str) -> float:
    return SLOW_INTERVAL if any(host == h or host.endswith("." + h) for h in SLOW_HOSTS) else HOST_INTERVAL


def _limit_file(host: str, suffix: str) -> Path:
    d = cache_dir() / "ratelimit"
    d.mkdir(exist_ok=True)
    return d / (re.sub(r"[^\w.-]", "_", host) + suffix)


def throttle(url: str, pause: float = 0.0) -> None:
    """Wait until this host's interval since the previous request (any process) has passed."""
    host = host_of(url)
    with open(_limit_file(host, ".lock"), "a+") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX)
        fh.seek(0)
        try:
            last = float(fh.read().strip() or 0)
        except ValueError:
            last = 0.0
        wait = max(pause, last + host_interval(host) - time.time())
        if wait > 0:
            time.sleep(wait)
        fh.seek(0)
        fh.truncate()
        fh.write(f"{time.time():.3f}")


def cooling_down(url: str) -> str:
    """-> why this host is in cool-down after a challenge/429, or ''."""
    path = _limit_file(host_of(url), ".blocked")
    try:
        since, why = path.read_text().split("\t", 1)
    except (OSError, ValueError):
        return ""
    left = float(since) + COOLDOWN_SECONDS - time.time()
    if left <= 0:
        path.unlink(missing_ok=True)
        return ""
    return f"host cool-down after {why.strip()} ({left / 3600:.1f} h left; not requested)"


def mark_blocked(url: str, why: str) -> None:
    _limit_file(host_of(url), ".blocked").write_text(f"{time.time():.0f}\t{why}")


def _curl(args: list[str]) -> subprocess.CompletedProcess:
    """Run curl; on an HTTP/2 stream error retry the same request once over HTTP/1.1."""
    done = subprocess.run(args, capture_output=True, check=False)
    if done.returncode == CURL_HTTP2_STREAM:
        done = subprocess.run(args[:1] + ["--http1.1"] + args[1:], capture_output=True, check=False)
    return done


def get(url: str, pause: float = 0.4, headers: dict | None = None,
        max_time: int = 120) -> tuple[int, bytes]:
    """-> (HTTP status, body) for small API calls. curl: the system Python has no CA bundle."""
    if cooling_down(url):
        return 0, b""
    throttle(url, pause)
    cmd = ["curl", "-sSL", "--connect-timeout", "30", "--max-time", str(max_time),
           "--max-filesize", str(MAX_BYTES), "-A", UA, "-w", "\n%{http_code}"]
    for k, v in (headers or {}).items():
        cmd += ["-H", f"{k}: {v}"]
    done = _curl(cmd + [url])
    body, _, code = done.stdout.rpartition(b"\n")
    status = int(code) if code.isdigit() else 0
    if status == 429:
        mark_blocked(url, "HTTP 429")
    else:
        gate = gate_reason(body[:65536])
        if gate and "Preparing to download" not in gate:
            mark_blocked(url, gate)
    return status, body


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
    blocked = cooling_down(url)
    if blocked:
        return 0, blocked
    throttle(url, 0.4)
    done = _curl(["curl", "-sSL", "--connect-timeout", "30", "--max-time",
                  str(max_time), "--speed-limit", "1", "--speed-time", str(STALL_SECONDS),
                  "--max-filesize", str(max_bytes), "-A", UA,
                  "-o", str(dest), "-w", "%{http_code}", url])
    out = done.stdout.decode(errors="replace").strip()
    code = int(out[-3:]) if out[-3:].isdigit() else 0
    if code == 429:
        mark_blocked(url, "HTTP 429")
    problem = {0: "", 63: f"larger than {max_bytes // 1_000_000} MB cap (not fetched)",
               28: "timed out or stalled (download incomplete)", 18: "transfer closed early (truncated)",
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
    if ERROR_BODY.search(head[:8192]):
        return f"error response, not a document ({ERROR_BODY.search(head).group(0).decode()})"
    if ext in PDF_EXT:
        if head.lstrip()[:5] != b"%PDF-" and b"%PDF-" not in head[:1024]:
            return gate or ("HTML page, not a PDF" if looks_html(head) else "not a PDF (no %PDF header)")
        if b"%%EOF" not in tail:
            return "truncated PDF (no %%EOF trailer)"
        return ""
    if ext in ZIP_EXT:
        if ext in OLE_FOR_OOXML and head[:8] == OLE_MAGIC:
            return ""       # a legacy Office file under an OOXML name (see content_name)
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
        if head[:8] == OLE_MAGIC or head[:2] == b"PK":
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

def content_name(name: str, path: Path) -> str:
    """The name a validated file is stored under: an OOXML name (.docx/.xlsx/.pptx) whose
    bytes are a legacy OLE2 Office file gets the legacy suffix, so readers pick the right
    converter. Every other name is returned unchanged."""
    p = Path(name)
    legacy = OLE_FOR_OOXML.get(p.suffix.lower())
    if legacy:
        with path.open("rb") as fh:
            if fh.read(8) == OLE_MAGIC:
                return p.stem + legacy
    return name


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

    def held(self, name: str) -> Path | None:
        """The file holding `name`'s content: `name` itself, or the byte-identical file a
        `(duplicate) name` row points to (by md5). None when neither is on disk."""
        if self.have(name):
            return self.path / name
        dup = self.rows.get(f"(duplicate) {name}")
        twin = self.by_md5.get(dup[3]) if dup and dup[3] else None
        if twin and (self.path / twin).exists():
            return self.path / twin
        return None

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
    def _url_tag(url: str) -> str:
        return hashlib.sha256(url.encode()).hexdigest()[:12]

    @classmethod
    def _manual_keys(cls, name: str, url: str) -> tuple[str, str]:
        path = Path(name)
        tagged = path.stem + "__" + cls._url_tag(url) + path.suffix
        return f"(manual) {name}", f"(manual) {tagged}"

    def _clear_manual(self, name: str, url: str) -> None:
        for key in self._manual_keys(name, url):
            if key in self.rows and self.rows[key][1] == url:
                self.rows.pop(key)

    def note(self, key: str, url: str, what: str) -> None:
        self.rows[f"(note) {key}"] = [f"(note) {key}", url, "", "", what]

    @staticmethod
    def mirror_name_kept(name: str, twin: str, twin_what: str) -> bool:
        """A scg_lib_structs repository copy of a vendor/publisher file (`scg_<file>`; not the
        repository's own scg_docs_/scg_adapters_ content) is kept under its scg_ name even
        when a byte-identical upstream_* file is on disk. upstream_* is the drawn page's
        namespace, which the source-transcription guard excludes; without the scg_ copy
        the original would be invisible to it."""
        return (name.startswith("scg_") and not name.startswith(("scg_docs_", "scg_adapters_"))
                and twin.startswith("upstream_") and twin_what != "scg_lib_structs method page")

    def _record(self, name: str, url: str, what: str) -> bool:
        """Row for a validated file now at self.path/name. False if it duplicates one."""
        p = self.path / name
        digest = md5(p)
        twin = self.by_md5.get(digest)
        if twin and twin != name and (self.path / twin).exists():
            trow = self.rows.get(twin)
            if self.mirror_name_kept(name, twin, trow[4] if trow else ""):
                what += f"; byte-identical to {twin}, kept under the repository-mirror name"
            else:
                p.unlink()
                if trow is not None and not trow[1] and url:
                    # a hand-placed original whose source URL is now known by its bytes
                    trow[1] = url
                    trow[4] += f"; URL back-filled from a byte-identical download ({what})"
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

    # --- source identity
    def _assigned_elsewhere(self, name: str, url: str) -> list[str] | None:
        """The row registering `name` to a different source URL, if any."""
        row = self.rows.get(name) or self.rows.get(f"(duplicate) {name}")
        return row if row is not None and row[1] and row[1] != url else None

    def _mirror_key(self, name: str, url: str) -> str:
        return f"(duplicate) {name}__{self._url_tag(url)}"

    def _collision_name(self, name: str, url: str) -> str:
        return "source-collision-" + name + "-" + self._url_tag(url)

    def _known_mirror(self, name: str, url: str) -> Path | None:
        """The registered file a previous run found byte-identical at `url`."""
        row = self.rows.get(self._mirror_key(name, url))
        held = self.held(name)
        if row and row[1] == url and row[3] and held and md5(held) == row[3]:
            return held
        return None

    def _differs_recorded(self, name: str, url: str) -> bool:
        return any(r[1] == url and "content differs" in r[4] for k, r in self.rows.items()
                   if k.startswith("(manual) " + self._collision_name(name, url)))

    def _mirror_or_collision(self, name: str, url: str, what: str, prior: list[str],
                             candidate: Path | None, problem: str) -> Path | None:
        """`name` is registered to another URL (`prior`). A validated `candidate` from `url`
        that is byte-identical to the registered file is the same original at another
        address (a publisher mirror): the alternate URL is recorded and the registered file
        returned. Anything else keeps the identity guard: the registered file is neither
        reused nor overwritten, and a `(manual)` row says why."""
        registered = self.held(name)
        digest = md5(candidate) if candidate else ""
        if candidate:
            candidate.unlink(missing_ok=True)
        if registered and digest and digest == md5(registered):
            key = self._mirror_key(name, url)
            self.rows[key] = [key, url, "", digest, f"{what}; byte-identical copy of "
                              f"{registered.name}, registered from {prior[1]}"]
            self._clear_manual(self._collision_name(name, url), url)
            self._clear_manual(name, url)
            say(f"    mirror  {name}  = {registered.name}  <- {url}")
            self.write()
            return registered
        detail = (f"content differs from the registered file (md5 {digest})" if digest else
                  f"identity not checkable ({problem or 'registered file missing'})")
        self.manual(self._collision_name(name, url), url, what,
                    f"cache filename {name} is already assigned to {prior[1]}; {detail}; "
                    "use a distinct source filename; existing file was not reused or overwritten")
        return None

    def _saved_from(self, url: str, name: str) -> str:
        """Another on-disk file already saved from exactly this URL (a second request for the
        same source would only differ by download stamping)."""
        for k, r in self.rows.items():
            if k != name and not k.startswith("(") and r[1] == url and (self.path / k).exists():
                return k
        return ""

    # --- fetching
    def _fetch(self, url: str, part: Path, name: str, min_text: int, need: bytes | None,
               max_bytes: int, retries: int) -> str:
        """Download `url` to `part` and validate it as `name`. -> '' or the reason it failed
        (then `part` is gone)."""
        host = re.sub(r"^\w+://([^/]+).*", r"\1", url)
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
                return ""
            gated = "gate" in reason or "challenge" in reason or "reCAPTCHA" in reason
            if gated and "Preparing to download" not in reason:
                # PMC's per-article download gate is not a host-wide bot challenge
                mark_blocked(url, reason)
            if (attempt == retries or gated or code == 429 or code not in RETRY_STATUS
                    or "cap" in reason or "stalled" in reason or "cool-down" in reason):
                break
            time.sleep(5 * (attempt + 1) ** 2)
        part.unlink(missing_ok=True)
        if "gate" in reason or "challenge" in reason or "reCAPTCHA" in reason or "cool-down" in reason:
            self.gated_hosts.add(host)
        return reason

    def _store(self, part: Path, name: str, url: str, what: str) -> Path | None:
        """Move a validated download into place under its content name and record it."""
        real = content_name(name, part)
        self._clear_manual(name, url)
        if real != name:
            prior = self._assigned_elsewhere(real, url)
            if prior is not None:    # the legacy file is already here from another URL
                held = self._mirror_or_collision(real, url, what, prior, part, "")
                if held:
                    self.rows[f"(duplicate) {name}"] = [f"(duplicate) {name}", url, "", md5(held),
                                                        f"{what}; legacy OLE2 content, identical to {held.name}"]
                    self.write()
                return held
        dest = self.path / real
        part.replace(dest)
        if real != name:
            label = f"{what} (legacy OLE2 Office file served under the name {name})"
            kept = self._record(real, url, label)
            digest = self.rows[real][3] if kept else self.rows[f"(duplicate) {real}"][3]
            self.rows[f"(duplicate) {name}"] = [f"(duplicate) {name}", url, "", digest,
                                                f"{what}; legacy OLE2 content stored as {real}"]
            self.write()
        else:
            kept = self._record(real, url, what)
        if not kept:
            return None
        say(f"    saved   {real}  ({dest.stat().st_size:,} bytes)")
        self.to_text(dest)
        return dest

    def save(self, name: str, url: str, what: str, *, min_text: int = 1500,
             need: bytes | None = None, max_bytes: int = MAX_BYTES,
             retries: int = 2) -> Path | None:
        """Download `url` as `name` unless already here; validate; convert to .txt.
        -> the saved path, or None (then a `(manual)` or `(duplicate)` row says why)."""
        name = safe_name(name)
        prior = self._assigned_elsewhere(name, url)
        if prior is not None:
            known = self._known_mirror(name, url)
            if known or self._differs_recorded(name, url):
                return known        # None: a previous run found different content there
            if self.held(name) is None:
                return self._mirror_or_collision(name, url, what, prior, None, "")
            part = self.path / (name + ".mirror.part")
            problem = self._fetch(url, part, name, min_text, need, max_bytes, retries=0)
            return self._mirror_or_collision(name, url, what, prior,
                                             None if problem else part, problem)
        dest = self.path / name
        if self.have(name) and not validate(dest, name, min_text, need):
            return dest
        if dest.exists() and not validate(dest, name, min_text, need):  # on disk, no row
            if self._record(name, url, what + " (already on disk)"):
                self.to_text(dest)
                return dest
            return None
        dup = self.rows.get(f"(duplicate) {name}")
        if dup:
            twin = self.by_md5.get(dup[3], "") if dup[3] else ""
            trow = self.rows.get(twin)
            if not (twin and self.mirror_name_kept(name, twin, trow[4] if trow else "")):
                if dup[1] == url and twin and (self.path / twin).exists():
                    self._clear_manual(name, url)     # an earlier failure at this URL is moot
                    self.write()
                return None
            del self.rows[f"(duplicate) {name}"]     # fetch again under the mirror name
        same = self._saved_from(url, name)
        if same:
            self.rows[f"(duplicate) {name}"] = [f"(duplicate) {name}", url, "", self.rows[same][3],
                                                f"{what}; same source URL already saved as {same}"]
            self._clear_manual(name, url)
            say(f"    dup     {name}  = {same} (same URL)")
            self.write()
            return None
        part = self.path / (name + ".part")
        reason = self._fetch(url, part, name, min_text, need, max_bytes, retries)
        if reason:
            self.manual(name, url, what, reason)
            return None
        dest = self._store(part, name, url, what)
        if dest and dest.suffix.lower() == ".zip" and self.expand:
            self.expand_zip(dest.name, dest.stem + "_", url, what)
        return dest

    def keep(self, name: str, data: bytes, url: str, what: str) -> Path | None:
        """Store bytes we already hold (from a zip, a cache) -- validated like a download."""
        name = safe_name(name)
        part = self.path / (name + ".part")
        part.write_bytes(data)
        reason = validate(part, name)
        prior = self._assigned_elsewhere(name, url)
        if prior is not None:
            if reason:
                part.unlink(missing_ok=True)
            return self._mirror_or_collision(name, url, what, prior,
                                             None if reason else part, reason)
        if reason:
            part.unlink()
            self.manual(name, url, what, reason)
            return None
        return self._store(part, name, url, what)

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
                        self._clear_manual(target, f"{url}#{m}")    # an earlier failure is moot
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
