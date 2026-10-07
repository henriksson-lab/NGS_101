#!/usr/bin/env python3
"""Render every .md note in the repo to a sibling .html, plus a notes index.

The protocol pages (atrandi_wgs.html, crisprmip.html, ...) are hand-built by each
protocol's tools/build_page.py because they carry sequence diagrams. Everything
else in the repo is Markdown, and this renders all of it through one shared
stylesheet so the notes and the diagrams look like the same publication.

    python3 build_docs.py            # render everything that changed
    python3 build_docs.py --force     # re-render regardless of timestamps
    python3 build_docs.py --check     # exit 1 if any .html is stale (for CI)
    python3 build_docs.py --list      # show what would be rendered

Markdown links to .md files are rewritten to .html, so the generated pages
navigate exactly like the sources.
"""
from __future__ import annotations

import argparse
import html
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "lib"))

from mdfacts import expand  # noqa: E402
from mdrender import render  # noqa: E402
from page import MD_STYLE, STYLE  # noqa: E402

# Directories that never contain notes to publish.
SKIP_DIRS = {".git", "_data", "pdf", "__pycache__", ".venv", "venv",
             "node_modules", ".ruff_cache", ".pytest_cache"}

# This one documents the upstream HTML conventions by showing raw markup; it is
# reference material for writing the diagram pages, not a note about chemistry.
INDEX_EXCLUDE = {"ref/scg_lib_structs_style.md"}

OUT_INDEX = ROOT / "notes.html"


def notes() -> list[Path]:
    found = [p for p in ROOT.rglob("*.md")
             if not SKIP_DIRS & set(p.relative_to(ROOT).parts)]
    return sorted(found, key=lambda p: str(p.relative_to(ROOT)).lower())


def title_of(md: str, path: Path) -> str:
    for line in md.split("\n"):
        if line.startswith("# "):
            # strip the inline markup a title might carry
            return html.unescape(plain(line[2:]))
    return path.stem.replace("_", " ").replace("-", " ")


def plain(text: str) -> str:
    """Heading text with inline markup removed, for a TOC entry or a title."""
    t = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", text)       # links -> their label
    t = re.sub(r"`|\*\*|\*|_", "", t)
    return html.escape(t.strip())


def toc_html(heads: list[tuple[int, str, str]]) -> str:
    """A table of contents, for notes long enough to need one."""
    items = [(lvl, txt, sid) for lvl, txt, sid in heads if lvl in (2, 3)]
    if len(items) < 4:
        return ""
    lis = []
    for lvl, txt, sid in items:
        # The TOC entry is itself a link, so the label must be plain text: a
        # heading containing a link would otherwise nest <a> inside <a>.
        lis.append(f'<li class="l{lvl}"><a href="#{sid}">{plain(txt)}</a></li>')
    return ('<nav class="toc">\n<b>On this page</b>\n<ul>\n'
            + "\n".join(lis) + "\n</ul>\n</nav>")


def build_page(path: Path) -> str:
    rel = path.relative_to(ROOT)
    md = expand(path.read_text(encoding="utf-8"), str(rel))
    depth = len(rel.parts) - 1
    up = "../" * depth
    heads: list[tuple[int, str, str]] = []
    body = render(md, headings=heads)
    title = title_of(md, path)

    crumb = (f'<nav class="crumb"><a href="{up}index.html">chem</a> &nbsp;/&nbsp; '
             f"{html.escape(str(rel))} &nbsp;&middot;&nbsp; "
             f'<a href="{up}notes.html">all notes</a> &nbsp;&middot;&nbsp; '
             f'<a href="{html.escape(path.name)}">Markdown source</a></nav>')
    foot = (f'<footer class="docfoot">Generated from <code>{html.escape(str(rel))}</code> '
            f"by <code>build_docs.py</code>. Edit the Markdown, not this file, then run "
            f"<code>python3 build_docs.py</code>.</footer>")
    return (f"<title>{html.escape(title)}</title>\n{STYLE}\n{MD_STYLE}\n"
            f'<div class="wrap">\n{crumb}\n{toc_html(heads)}\n'
            f'<article class="md">\n{body}\n</article>\n{foot}\n</div>\n')


def build_index(paths: list[Path]) -> str:
    groups: dict[str, list[Path]] = {}
    for p in paths:
        rel = p.relative_to(ROOT)
        if str(rel).replace("\\", "/") in INDEX_EXCLUDE:
            continue
        groups.setdefault(rel.parts[0] if len(rel.parts) > 1 else ".", []).append(p)

    secs = []
    for key in sorted(groups, key=lambda k: (k == ".", k)):
        rows = []
        for p in groups[key]:
            rel = p.relative_to(ROOT)
            md = p.read_text(encoding="utf-8")
            lines = len(md.split("\n"))
            rows.append(
                f'<a href="{html.escape(str(rel.with_suffix(".html")))}">'
                f'<span class="n">{html.escape(title_of(md, p))}</span>'
                f'<span class="m">{lines:,} lines</span></a>')
        name = "repository root" if key == "." else key
        secs.append(f"<h2>{html.escape(name)}</h2>\n"
                    f'<div class="doclist">\n' + "\n".join(rows) + "\n</div>")

    total = sum(len(v) for v in groups.values())
    return ("<title>Notes index</title>\n" + STYLE + "\n" + MD_STYLE + "\n"
            '<div class="wrap">\n'
            '<nav class="crumb"><a href="index.html">chem</a> &nbsp;/&nbsp; notes</nav>\n'
            "<h1>Notes</h1>\n"
            f'<article class="md">\n<p>Every Markdown note in the repository, rendered '
            f"by <code>build_docs.py</code> &mdash; {total} notes. The protocol pages with "
            f'sequence diagrams are linked from the <a href="index.html">front page</a> '
            f"instead, since those are built separately.</p>\n"
            + "\n".join(secs) + "\n</article>\n</div>\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true",
                    help="do not write; exit 1 if any output is missing or stale")
    ap.add_argument("--force", action="store_true", help="re-render everything")
    ap.add_argument("--list", action="store_true", help="list the notes and exit")
    ap.add_argument("paths", nargs="*", type=Path, help="limit to these .md files")
    args = ap.parse_args()

    found = notes()
    if args.paths:
        want = {p.resolve() for p in args.paths}
        found = [p for p in found if p.resolve() in want]
        if not found:
            print("no matching .md files", file=sys.stderr)
            return 2

    if args.list:
        for p in found:
            print(p.relative_to(ROOT))
        return 0

    targets = [(p, p.with_suffix(".html")) for p in found]
    targets.append((None, OUT_INDEX))
    stale, written = [], 0

    for src, dst in targets:
        new = build_index(found) if src is None else build_page(src)
        old = dst.read_text(encoding="utf-8") if dst.exists() else None
        if old == new and not args.force:
            continue
        if args.check:
            stale.append(dst.relative_to(ROOT))
            continue
        dst.write_text(new, encoding="utf-8")
        written += 1
        print(f"  {dst.relative_to(ROOT)}")

    if args.check:
        if stale:
            print(f"{len(stale)} page(s) stale -- run: python3 build_docs.py",
                  file=sys.stderr)
            for s in stale:
                print(f"  {s}", file=sys.stderr)
            return 1
        print(f"all {len(targets)} page(s) up to date")
        return 0

    print(f"{written} written, {len(targets) - written} unchanged "
          f"({len(found)} notes + notes.html)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
