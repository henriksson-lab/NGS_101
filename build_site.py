#!/usr/bin/env python3
"""Build the website: every page, then a deployable `_site/` holding only what we wrote.

    python3 build_site.py                 # build everything, assemble and check _site/
    python3 build_site.py --no-checks     # skip running the self-tests for the index
    python3 -m http.server -d _site       # preview at http://localhost:8000/

Steps:

1. Run every `*/tools/build_page.py`. A failure is reported, not fatal: in CI the
   third-party source files (plasmid maps, supplementary tables) are absent, and those
   build scripts exit with a "missing source" message. Such a page is left out of the site
   -- the index does not link it, and a short stand-in page explains why, so that a link
   to it from another page still lands somewhere.
2. Run `build_docs.py` (maintainer-facing rendered notes) and `build_index.py`.
3. Assemble `_site/` from a whitelist: index.html, the diagram pages that built, each
   diagram's main research note, and files those pages link to -- provided they are ours.
   Notes therefore appear only as supporting material reached from a finished schematic.
   Never third-party material:
   nothing from `_data/`, `pdf/`, download
   caches, archived exemplars (`ref/*.html` with no Markdown twin), or any `ref/`
   directory's data files, and no file of a type `.gitignore` treats as source material.
4. Check every relative link and anchor in `_site/`; a broken one fails the build.
"""
from __future__ import annotations

import argparse
import html
import os
import posixpath
import re
import shutil
import subprocess
import sys
import urllib.parse
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SITE = ROOT / "_site"
sys.path.insert(0, str(ROOT / "lib"))
sys.path.insert(0, str(ROOT / "catalogue" / "tools"))
sys.path.insert(0, str(ROOT))

import build_docs  # noqa: E402
import build_index  # noqa: E402
from page import ASSETS, MD_STYLE  # noqa: E402

# Never published from here, whatever links to it.
NEVER = {"_data", "pdf", ".cache", ".git", "__pycache__", "_site"}
# Files a page may link to and that we may copy, if they are ours. Everything else -- in
# particular the types .gitignore lists as third-party source material (pdf, xlsx, txt,
# images, plasmid maps) -- is never copied.
ASSET_TYPES = {".py", ".R", ".sh", ".css", ".js", ".svg", ".json", ".tsv", ".csv"}
THIRD_PARTY_TYPES = {".pdf", ".xls", ".xlsx", ".doc", ".docx", ".pptx", ".txt", ".jpg",
                     ".jpeg", ".png", ".gif", ".tif", ".tiff", ".dna", ".gb", ".gbk",
                     ".fa", ".fasta"}


def rel(p: Path) -> str:
    return p.relative_to(ROOT).as_posix()


# ----------------------------------------------------------------- 1. build
def build_pages() -> dict[str, tuple[bool, str]]:
    """Run each protocol's build_page.py; {dir: (ok, last line of its output)}."""
    out = {}
    for bp in sorted(ROOT.glob("*/tools/build_page.py")):
        d = bp.parents[1].name
        r = subprocess.run([sys.executable, str(bp)], capture_output=True, text=True,
                           cwd=ROOT)
        msg = (r.stdout + r.stderr).strip().splitlines()
        last = msg[-1] if msg else ""
        out[d] = (r.returncode == 0, last)
        print(f"  {'ok  ' if r.returncode == 0 else 'FAIL'}  {d}: {last[:150]}")
    return out


def run(script: str, *args: str) -> None:
    r = subprocess.run([sys.executable, str(ROOT / script), *args], cwd=ROOT)
    if r.returncode:
        sys.exit(f"build_site: {script} failed (exit {r.returncode})")


# -------------------------------------------------------------- 3. assemble
def published(p: Path) -> bool:
    """May this repository file appear in the site at all?"""
    parts = p.relative_to(ROOT).parts
    if NEVER & set(parts) or p.suffix.lower() in THIRD_PARTY_TYPES:
        return False
    if "ref" in parts and p.suffix not in (".md", ".html"):
        return False                  # ref/ data files are archived third-party material
    if p.suffix == ".html" and not p.with_suffix(".md").exists():
        return False                  # only rendered notes; diagram pages are added by name
    return True


def placeholder(d: str, out: Path, why: str) -> str:
    up = "../" * (len(out.relative_to(ROOT).parts) - 1)
    return (f"<!doctype html>\n<title>{html.escape(out.name)} — not in this build</title>\n"
            f"{ASSETS}\n{MD_STYLE}\n<div class=\"wrap\">\n"
            f'<nav class="crumb"><a href="{up}index.html">chem</a> &nbsp;/&nbsp; '
            f"{html.escape(d)}</nav>\n<article class=\"md\">\n"
            f"<h1>{html.escape(out.name)} is not part of this build</h1>\n"
            f"<p>This diagram page is generated from third-party source files (papers, "
            f"supplementary tables, plasmid maps) that are not redistributed with the "
            f"repository, so the public build cannot draw it. Build it locally once the "
            f"files listed in the protocol's <code>ref/MANIFEST.md</code> are in place.</p>\n"
            f"<p>The build script said: <code>{html.escape(why)}</code></p>\n"
            + "</article>\n</div>\n")


class Links(HTMLParser):
    """Every href/src on a page, plus every id it defines."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.refs: list[str] = []
        self.ids: set[str] = set()

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if a.get("id"):
            self.ids.add(a["id"])
        if tag == "a" and a.get("name"):
            self.ids.add(a["name"])
        for k in ("href", "src"):
            if a.get(k) and tag in ("a", "link", "script", "img", "iframe", "source"):
                self.refs.append(a[k])


def parse(p: Path) -> Links:
    lp = Links()
    lp.feed(p.read_text(encoding="utf-8"))
    return lp


def local(ref: str) -> tuple[str, str] | None:
    """(path, fragment) of a relative link; None for external or scheme links."""
    u = urllib.parse.urlsplit(ref)
    if u.scheme or u.netloc or ref.startswith("//"):
        return None
    return urllib.parse.unquote(u.path), urllib.parse.unquote(u.fragment)


def resolve(page: str, path: str) -> str:
    """Site-relative target of `path` linked from site-relative `page`."""
    if not path:
        return page
    return posixpath.normpath(posixpath.join(posixpath.dirname(page), path))


H1_END = re.compile(r"(</h1>)", re.I)


def link_research_note(page: str, note: Path) -> str:
    """Add one quiet route from a schematic to its protocol's main research note."""
    if 'class="research-notes"' in page:
        return page
    link = (f'\n<p class="research-notes"><a href="{html.escape(note.with_suffix(".html").name)}">'
            'Research notes</a></p>')
    linked, n = H1_END.subn(r"\1" + link, page, count=1)
    if n != 1:
        raise ValueError(f"cannot add research-note link: schematic has {n} closing h1 tags")
    return linked


def assemble(results: dict[str, tuple[bool, str]]) -> list[str]:
    if SITE.exists():
        shutil.rmtree(SITE)
    SITE.mkdir()
    (SITE / ".nojekyll").write_text("", encoding="utf-8")
    files: dict[str, Path | str] = {}            # site path -> source file, or content

    files["index.html"] = ROOT / "index.html"

    skipped = []
    public_dirs = {r["dir"] for r in build_index.cat.ours() if r["section"] == "published"}
    for d, (ok, why) in sorted(results.items()):
        if d not in public_dirs:
            continue
        out = build_index.diagram_out(d)
        if out is None:
            continue
        if ok and out.exists():
            page = out.read_text(encoding="utf-8")
            notes = sorted((ROOT / d).glob("0*.md"))
            if notes:
                note = notes[0]
                note_html = note.with_suffix(".html")
                files[rel(note_html)] = note_html
                page = link_research_note(page, note)
            files[rel(out)] = page
        else:
            skipped.append(d)
            files[rel(out)] = placeholder(d, out, why or "build failed")

    # the files those pages link to, if they are ours and may be published
    pending = [k for k in files if k.endswith(".html")]
    seen = set(pending)
    while pending:
        page = pending.pop()
        src = files[page]
        text = src if isinstance(src, str) else src.read_text(encoding="utf-8")
        lp = Links()
        lp.feed(text)
        for ref in lp.refs:
            loc = local(ref)
            if not loc or not loc[0]:
                continue
            target = resolve(page, loc[0])
            if target in files or target.startswith("../"):
                continue
            p = ROOT / target
            if p.is_file() and p.suffix in ASSET_TYPES and published(p):
                files[target] = p
            elif p.is_dir() and (p / "README.md").exists() and published(p / "README.md"):
                # a link to a directory: serve its rendered README as the directory index
                files[target.rstrip("/") + "/index.html"] = p / "README.html"

    for k, src in files.items():
        dst = SITE / k
        dst.parent.mkdir(parents=True, exist_ok=True)
        if k.endswith(".html"):
            text = src if isinstance(src, str) else src.read_text(encoding="utf-8")
            text = MD_SOURCE.sub("", text)
            dst.write_text(unlink_unpublished(k, text, files), encoding="utf-8")
        else:
            shutil.copyfile(src, dst)
    return skipped


# A note page's breadcrumb links its Markdown source, which is useful in a checkout but
# not on the website: the site carries rendered pages only.
MD_SOURCE = re.compile(r' &nbsp;&middot;&nbsp; <a href="[^"/]*\.md">Markdown source</a>')
A_TAG = re.compile(r'<a href="([^"]*)"([^>]*)>(.*?)</a>', re.S)


def unlink_unpublished(page: str, text: str, files: dict) -> str:
    """Turn a link to something that exists in the repository but is deliberately not
    published (a debug log, a ref/ data file, a bare directory) into plain text with a
    tooltip. A link to something that does not exist at all is left alone, so that the
    link check still catches it."""
    def fix(m: re.Match) -> str:
        loc = local(html.unescape(m.group(1)))
        if not loc or not loc[0]:
            return m.group(0)
        target = resolve(page, loc[0])
        if target in files or target.rstrip("/") + "/index.html" in files:
            return m.group(0)
        if target.startswith("../") or not (ROOT / target).exists():
            return m.group(0)
        return (f'<span class="unpub" title="in the repository, not published on the '
                f'website: {html.escape(target)}" style="border-bottom:1px dotted '
                f'var(--ink-muted)">{m.group(3)}</span>')
    return A_TAG.sub(fix, text)


# ------------------------------------------------------------------ 4. check
def check_links() -> list[str]:
    """Every relative link in _site/ must resolve to a file in _site/ (and anchor)."""
    pages = {p.relative_to(SITE).as_posix(): p for p in SITE.rglob("*.html")}
    parsed = {k: parse(p) for k, p in pages.items()}
    bad = []
    for page, lp in sorted(parsed.items()):
        for ref in lp.refs:
            loc = local(ref)
            if loc is None:
                continue
            path, frag = loc
            target = resolve(page, path)
            if target.startswith("../") or target == "..":
                bad.append(f"{page}: {ref} -> leaves the site")
                continue
            tp = SITE / target
            if tp.is_dir():
                target = posixpath.join(target, "index.html")
                tp = SITE / target
            if not tp.is_file():
                bad.append(f"{page}: {ref} -> {target} does not exist")
            elif frag and target in parsed and frag not in parsed[target].ids:
                bad.append(f"{page}: {ref} -> no id {frag!r} in {target}")
    return bad


def audit() -> list[str]:
    """Nothing but our own output types made it into _site/."""
    problems = []
    for p in SITE.rglob("*"):
        if p.is_dir():
            continue
        r = p.relative_to(SITE)
        if p.suffix.lower() in THIRD_PARTY_TYPES or NEVER & set(r.parts):
            problems.append(f"{r}: third-party or unpublished material")
        if p.suffix not in ASSET_TYPES | {".html", ""}:
            problems.append(f"{r}: unexpected file type")
    return problems


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--no-checks", action="store_true",
                    help="do not run the self-tests to put check counts on the index")
    args = ap.parse_args()

    print("1. diagram pages (*/tools/build_page.py)")
    results = build_pages()
    failed = sorted(d for d, (ok, _) in results.items() if not ok)

    print("2. notes (build_docs.py) and index (build_index.py)")
    run("build_docs.py")
    run("build_index.py", *(a for d in failed for a in ("--omit", d)),
        *(["--no-checks"] if args.no_checks else []))

    print("3. assembling _site/")
    skipped = assemble(results)
    n = sum(1 for p in SITE.rglob("*") if p.is_file())
    size = sum(p.stat().st_size for p in SITE.rglob("*") if p.is_file())
    print(f"  {n} files, {size / 1e6:.1f} MB")

    print("4. checking links")
    bad = check_links() + audit()
    for b in bad:
        print(f"  BROKEN  {b}")

    built = len(results) - len(failed)
    public_dirs = {r["dir"] for r in build_index.cat.ours() if r["section"] == "published"}
    public_built = sum(d in public_dirs and ok for d, (ok, _) in results.items())
    print(f"\n{built} of {len(results)} diagram pages built; "
          f"{public_built} finished published schematics included"
          + (f"; left out (stand-in page instead): {', '.join(skipped)}" if skipped else ""))
    if bad:
        print(f"{len(bad)} problem(s) in _site/ -- not deployable", file=sys.stderr)
        return 1
    print(f"_site/ is ready: every relative link resolves. Preview: "
          f"python3 -m http.server -d {os.path.relpath(SITE)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
