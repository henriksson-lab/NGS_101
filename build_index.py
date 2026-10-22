#!/usr/bin/env python3
"""Generate index.html -- the public landing page of finished protocol schematics.

Nothing about an individual protocol is typed here. The list comes from
`catalogue/ours.tsv` (one row per directory; its `section` column says whether it is a
published protocol or our own work in progress), joined to the scraped catalogue for
category, family members, papers and year. Every published protocol is listed by name; its
name becomes a link only with a finished catalogue status and a generated diagram. Notes and work in progress
stay out of the public index.

The published schematics get a client-side property finder backed by the validated chemistry
facets in `catalogue/properties.tsv`; free text is secondary. No note bodies are embedded in
the public page. No library, no fetch() -- it works from file:// too, and without JavaScript
the full list is simply shown.

Check counts are read by actually running each suite, so the page cannot drift from
reality: if a suite is failing, the page says so.

    python3 build_index.py                  # write index.html
    python3 build_index.py --omit DIR ...   # treat DIR's diagram page as not built
    python3 build_index.py --no-checks      # skip running the self-tests (faster)
"""
from __future__ import annotations

import argparse
import html
import json
import re
import subprocess
import sys
from html.parser import HTMLParser
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "lib"))
sys.path.insert(0, str(ROOT / "catalogue" / "tools"))

import catalogue as cat  # noqa: E402
import citations as cites  # noqa: E402
import paper_overrides  # noqa: E402
import properties as props  # noqa: E402
import vendors  # noqa: E402

# File discovery also finds unfinished pages; it cannot promote catalogue status.
FINISHED_SCHEMATIC_STATUSES = frozenset({"documented", "assay", "construct"})
from mdfacts import expand  # noqa: E402
from mdrender import render  # noqa: E402
from page import head, info, site_header  # noqa: E402
from protocol_pages import page_module  # noqa: E402

OUT = ROOT / "index.html"

# Suites that are not a protocol directory but still count towards the total.
EXTRA_SUITES = ["catalogue"]

BLURB_CHARS = 330          # a blurb is cut at a sentence end before this many characters
# Build scripts all declare an OUT path, but older pages use either HERE.parent or
# Path(__file__).resolve().parents[1], with varying whitespace.  The filename after the
# final path-join slash is the stable part; page discovery must not depend on source style.
OUT_RE = re.compile(
    r'(?:^|;)\s*OUT\s*=[^\n]*/\s*(?P<quote>["\'])(?P<out>[^"\']+\.html)(?P=quote)',
    re.M,
)


# ------------------------------------------------------------------ the notes
class _Text(HTMLParser):
    """Rendered note -> (plain text of every block, top-level paragraphs in order)."""

    BLOCKS = {"p", "li", "td", "th", "h1", "h2", "h3", "h4", "h5", "h6", "pre",
              "blockquote", "tr", "div", "br"}

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.text: list[str] = []
        self.paras: list[tuple[str, str]] = []     # (h2 above it, paragraph text)
        self._quote = 0          # inside a blockquote
        self._quote_at = 0       # where in self.text the outermost blockquote began
        self._list = 0           # inside a list or table
        self._para: list[str] | None = None
        self._h2: list[str] | None = None
        self.section = ""        # text of the h2 the parser is under ("" before the first)

    def handle_starttag(self, tag, attrs):
        if tag == "blockquote":
            if not self._quote:
                self._quote_at = len(self.text)
            self._quote += 1
        elif tag in ("ul", "ol", "table"):
            self._list += 1
        elif tag == "p" and not self._quote and not self._list:
            self._para = []
        elif tag == "h2":
            self._h2 = []
        if tag in self.BLOCKS:
            self.text.append(" ")

    def handle_endtag(self, tag):
        if tag == "blockquote":
            self._quote -= 1
            # the evidence-marking key every note opens with is boilerplate, not content
            if not self._quote and _squash("".join(self.text[self._quote_at:])).startswith(
                    "Evidence marking"):
                del self.text[self._quote_at:]
        elif tag in ("ul", "ol", "table"):
            self._list -= 1
        elif tag == "p" and self._para is not None:
            self.paras.append((self.section, _squash("".join(self._para))))
            self._para = None
        elif tag == "h2" and self._h2 is not None:
            self.section = _squash("".join(self._h2))
            self._h2 = None
        if tag in self.BLOCKS:
            self.text.append(" ")

    def handle_data(self, data):
        self.text.append(data)
        if self._para is not None:
            self._para.append(data)
        if self._h2 is not None:
            self._h2.append(data)


def _squash(s: str) -> str:
    return re.sub(r"\s+", " ", s).strip()


def read_note(path: Path) -> dict:
    """Title, blurb and plain text of one Markdown note, facts expanded."""
    rel = path.relative_to(ROOT)
    md = expand(path.read_text(encoding="utf-8"), str(rel))
    title = ""
    for line in md.split("\n"):
        if line.startswith("# "):
            title = _squash(re.sub(r"\[([^\]]*)\]\([^)]*\)|[`*_]", lambda m: m.group(1) or "",
                                   line[2:]))
            break
    p = _Text()
    p.feed(render(md))
    return {"title": title or path.stem.replace("_", " "),
            "href": str(rel.with_suffix(".html")).replace("\\", "/"),
            "blurb": blurb(p.paras),
            "text": _squash("".join(p.text))}


# Paragraphs about the sources rather than the method: what was fetched, what could not be,
# which other catalogue rows exist. A blurb should say what the method IS.
NOT_BLURB = re.compile(
    r"\bdoi\b|fetch|obtained|download|paywall|transcribed from|catalogue|manifest|_data|"
    r"get_sources|^(sources?|status|other|related|original method|note)\b", re.I)
MARKS = re.compile(r"\s*[\U0001F7E2\U0001F7E1\U0001F534\u2705\u26A0]\uFE0F?")


# The section most notes open with ("## 1. What it is"), whose first paragraph is the
# best one-paragraph summary of the method; failing that, the note's own introduction.
# ...and paragraphs that only make sense after the one before them.
CONTINUES = re.compile(r"(so|thus|hence|this note|this one note|both are)\b", re.I)
SUMMARY_H2 = re.compile(r"what it is|what it replaces|what changed|principle|conceptual point",
                        re.I)


def blurb(paras: list[tuple[str, str]]) -> str:
    """The first paragraph that says what the thing IS -- not a citation or a file list."""
    def ok(t: str) -> bool:
        return not (len(t) < 60 or NOT_BLURB.search(t) or t.endswith((":", "\u2014"))
                    or t[0] in "(" or CONTINUES.match(t)
                    or "therefore" in t.split(". ")[0])

    def lead_in(t: str) -> str:
        """A paragraph introducing a list still reads as a summary, ended with a stop."""
        return t[:-1].rstrip() + "." if t.endswith(":") else t

    cands = [lead_in(MARKS.sub("", t).strip()) for _, t in paras]
    order = ([c for (h, _), c in zip(paras, cands) if SUMMARY_H2.search(h)]
             + [c for (h, _), c in zip(paras, cands) if not h] + cands)
    for t in order:
        if not t or not ok(t):
            continue
        if len(t) <= BLURB_CHARS:
            return t
        cut = t[:BLURB_CHARS]
        end = max(cut.rfind(". "), cut.rfind("; "))
        return (cut[:end + 1] if end > 80 else cut.rsplit(" ", 1)[0]) + " \u2026"
    return ""


def notes_of(d: str) -> list[Path]:
    return sorted((ROOT / d).glob("0*.md"))


# ---------------------------------------------------------- the diagram pages
def diagram_out(d: str) -> Path | None:
    """Output of a legacy page builder or a PAGE_FILE protocol module."""
    bp = ROOT / d / "tools" / "build_page.py"
    if bp.exists():
        m = OUT_RE.search(bp.read_text(encoding="utf-8"))
        if not m:
            raise ValueError(f"{bp}: cannot determine generated HTML from its OUT assignment")
        return ROOT / d / m.group("out")
    declared = page_module(ROOT / d)
    return ROOT / d / declared[1] if declared else None


def checks(d: str) -> tuple[int, bool]:
    """Run a directory's self-test and report (count, passed); (0, False) if none."""
    st = ROOT / d / "tools" / "selftest.py"
    if st.exists():
        command = [sys.executable, str(st)]
    elif page_module(ROOT / d):
        command = [sys.executable, str(ROOT / "tools" / "protocol.py"), "test", d]
    else:
        return 0, False
    r = subprocess.run(command, capture_output=True, text=True)
    m = re.search(r"all (\d+) checks passed", r.stdout)
    if m and r.returncode == 0:
        return int(m.group(1)), True
    m = re.search(r"(\d+) of (\d+) checks FAILED", r.stdout)
    return (int(m.group(2)) if m else 0), False


# ------------------------------------------------------------- one protocol
CATEGORY_LABEL = {"TODO list": "on the scg_lib_structs TODO list",
                  "Ours, not in scg_lib_structs": "not in scg_lib_structs"}


def entry(r: dict, properties: dict, citation: cites.Citation | None,
          omit: set[str], run_checks: bool) -> dict:
    d = r["dir"]
    crows = cat.rows_for_dir(d)
    # Our one-protocol-per-paper record is authoritative when an upstream aggregate
    # (for example the former SMART-seq family page) marks a different family paper as
    # its defining row.
    defining = (next((x for x in crows if x["doi"].lower() == r["doi"].lower()
                      and r["doi"]), None)
                or next((x for x in crows if x["is_defining"] == "yes"), {}))
    primary = [x for x in crows if x["role"] == "primary"]
    # Other names it goes by: family members, and the catalogue's protocol name where it
    # differs from ours. (A paper's anchor text is searchable but not shown: upstream
    # anchors are as often "Nature Protocols 12, 44-73" as a method name.)
    aka: list[str] = []
    for x in crows:
        for name in [x["protocol"], *x["family_members"].split(";")]:
            name = name.strip()
            if name and name.lower() != r["protocol"].lower() \
                    and name.lower() not in {a.lower() for a in aka}:
                aka.append(name)
    notes = [read_note(p) for p in notes_of(d)]
    page = diagram_out(d)
    has_page = bool(r["status"] in FINISHED_SCHEMATIC_STATUSES
                    and page and page.exists() and d not in omit)
    if has_page and r["status"] in {"documented", "construct"}:
        has_page = 'data-protocol-complete="true"' in page.read_text(encoding="utf-8")
        if has_page and r["status"] == "construct":
            has_page = 'data-protocol-output="constructed DNA"' in page.read_text(encoding="utf-8")
    n, ok = checks(d) if run_checks else (0, False)
    category = crows[0]["category"] if crows else ""
    commercial = "commercial kit" in properties["availability"]
    if r["doi"]:
        override = paper_overrides.BY_DIRECTORY.get(d)
        author = (override.first_author if override else
                  html.unescape(defining.get("first_author", "")).strip())
        year = override.year if override else defining.get("year", "").strip()
        journal = (override.journal if override else
                   html.unescape(defining.get("journal", "")).strip())
        if not (author and year and journal):
            raise ValueError(f"{d}: defining paper lacks author/year/journal metadata")
        reference = f"{author} et al., {year}, {journal}"
        reference_href = override.href if override else f"https://doi.org/{r['doi']}"
    elif commercial:
        reference = vendors.company_for(d)
        reference_href = ""
    else:
        reference = ""
        reference_href = ""
    return {
        "dir": d, "name": r["protocol"], "doi": r["doi"], "status": r["status"],
        "section": r["section"], "modality": r["modality"],
        "category": CATEGORY_LABEL.get(category, category),
        "aka": aka,
        "year": defining.get("year", ""),
        "papers": sorted({x["title"] for x in primary if x["title"]}
                         | {x["citation_text"] for x in primary if x["citation_text"]}),
        "scg_page": next((x["scg_page"] for x in crows if x["scg_page"]), ""),
        "notes": [{"title": n_["title"], "href": n_["href"]} for n_ in notes],
        "blurb": next((n_["blurb"] for n_ in notes if n_["blurb"]), r["note"]),
        "text": " ".join(n_["text"] for n_ in notes),
        "page": (str(page.relative_to(ROOT)).replace("\\", "/") if has_page else ""),
        "checks": n, "checks_ok": ok,
        "properties": properties,
        "citations": citation.citations if citation else None,
        "citations_retrieved": citation.retrieved if citation else "",
        "citation_openalex": citation.openalex_id if citation else "",
        "commercial": commercial,
        "reference": reference,
        "reference_href": reference_href,
    }


def collect(omit: set[str] = frozenset(), run_checks: bool = True) -> list[dict]:
    ours = cat.ours()
    properties = props.load({r["dir"] for r in ours})
    citations = cites.load(expected_dois={r["doi"].strip().lower()
                                          for r in ours if r["doi"].strip()})
    return [entry(r, properties[r["dir"]], citations.get(r["doi"].lower()),
                  set(omit), run_checks) for r in ours]


# ------------------------------------------------------------------- render
e = html.escape


def status_label(p: dict) -> str:
    return "schematic" if p["page"] else "source notes"


def links(p: dict) -> str:
    out = []
    if p["page"]:
        out.append(f'<a class="go" href="{e(p["page"])}">open schematic &rarr;</a>')
    if p["scg_page"]:
        out.append(f'<a class="ext source" href="{e(p["scg_page"])}">original schematic</a>')
    return '<div class="links">' + "".join(out) + "</div>"


def meta_line(p: dict) -> str:
    bits = []
    if p["year"]:
        bits.append(e(p["year"]))
    if p["doi"]:
        bits.append(f'<a href="https://doi.org/{e(p["doi"])}">doi:{e(p["doi"])}</a>')
    return '<div class="meta">' + " &middot; ".join(bits) + "</div>"


def badges(p: dict) -> str:
    b = [f'<span class="tag mod">{e(p["modality"])}</span>']
    if p["category"]:
        b.append(f'<span class="tag">{e(p["category"])}</span>')
    return '<div class="tags">' + "".join(b) + "</div>"


def result_item(i: int, p: dict) -> str:
    first = p["page"]
    name = f'<a href="{e(first)}">{e(p["name"])}</a>' if first else e(p["name"])
    if p["citations"] is not None:
        n = p["citations"]
        citation = (f'<a class="cites" href="{e(p["citation_openalex"])}" '
                    f'title="OpenAlex snapshot, {e(p["citations_retrieved"])}">'
                    f'{n:,}</a>')
    elif not p["doi"] and p["commercial"]:
        citation = '<span class="cites no-count">commercial protocol</span>'
    elif not p["doi"]:
        citation = '<span class="cites no-count">no defining paper</span>'
    else:
        citation = '<span class="cites no-count">citation count unavailable</span>'
    family = ", ".join(p["properties"]["assay"])
    scale = ", ".join(p["properties"]["scale"])
    platform = ", ".join(p["properties"]["platform"])
    if p["reference_href"]:
        reference = f'<a href="{e(p["reference_href"])}">{e(p["reference"])}</a>'
    elif p["reference"]:
        reference = f'<span class="company">{e(p["reference"])}</span>'
    else:
        reference = '<span class="no-reference">&mdash;</span>'
    return f"""<tr class="pr" data-i="{i}" id="p-{e(p['dir'])}">
<th scope="row">{name}</th>
<td class="reference-cell">{reference}</td>
<td class="citation-cell">{citation}</td>
<td>{e(family)}</td>
<td>{e(scale)}</td>
<td>{e(platform)}</td>
</tr>"""


def search_data(ps: list[dict]) -> str:
    """The search index: one record per published protocol, in list order."""
    recs = [{"n": p["name"], "a": p["aka"], "c": p["category"], "m": p["modality"],
             "s": status_label(p),
             "d": p["doi"], "y": p["year"],
             "r": p["reference"],
             "t": " ".join([*p["papers"], p["blurb"]]),
             "p": p["properties"], "x": p["citations"]} for p in ps]
    # embedded in a <script>, so no "</" may appear literally
    return json.dumps(recs, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")


EXTRA_CSS = """<style>
.sr-only { position:absolute; width:1px; height:1px; padding:0; margin:-1px;
           overflow:hidden; clip:rect(0,0,0,0); white-space:nowrap; border:0; }
.cards { display:grid; gap:14px; grid-template-columns:repeat(auto-fit,minmax(min(330px,100%),1fr));
         margin:1.2em 0; align-items:start; }
.card { border:1px solid var(--rule); border-radius:6px; padding:16px 18px;
        background:var(--surface); display:flex; flex-direction:column; gap:8px; }
.card h3, .concept h3 { margin:0; font-size:1.06rem; }
.card h3::before, .concept h3::before { content:none; }
.card p { margin:0; font-size:.95rem; }
.kind { font-size:.72rem; letter-spacing:.08em; text-transform:uppercase;
        color:var(--ink-muted); font-weight:600; }
.tags { display:flex; flex-wrap:wrap; gap:5px; }
.tag { font-size:.76rem; padding:2px 8px; border-radius:999px;
       background:var(--surface-2); color:var(--ink-muted); }
.tag.mod { font-weight:600; color:var(--ink); }
.tag.pg { color:var(--c-tso); }
.meta { font-size:.8rem; color:var(--ink-muted); font-variant-numeric:tabular-nums;
        overflow-wrap:anywhere; }
.links { display:flex; flex-wrap:wrap; gap:4px 14px; font-size:.9rem; }
.links .go { font-weight:600; }
.links .ext::after { content:' \\2197'; }
.wipbox { border-left:3px solid var(--note-rule); padding-left:14px; }

/* property finder */
.index { padding-top:14px; }
.index .wrap > h1 { font-size:1.55rem; margin:0; line-height:1.1; }
.finder { margin:.35em 0 .55em; padding:8px 10px; border:1px solid var(--rule); border-radius:6px;
          background:var(--surface); }
.finder-head { display:flex; gap:8px; align-items:baseline; margin-bottom:5px; }
.finder-head h2 { margin:0; font-size:1rem; }
.finder-head h2::before { content:none; }
.facetbar { display:flex; flex-wrap:wrap; gap:5px; align-items:flex-start; }
.facet { position:relative; }
.facet summary, .more-facets > summary, .text-search > summary { list-style:none; cursor:pointer;
          font-size:.8rem; padding:3px 8px; border:1px solid var(--rule); border-radius:999px;
          color:var(--ink-muted); background:var(--bg); user-select:none; }
.facet summary::-webkit-details-marker, .more-facets > summary::-webkit-details-marker,
.text-search > summary::-webkit-details-marker { display:none; }
.facet summary::after, .more-facets > summary::after, .text-search > summary::after {
          content:' +'; color:var(--accent); }
.facet[open] summary, .facet summary.active { border-color:var(--accent); color:var(--ink); }
.facet[open] summary::after, .more-facets[open] > summary::after,
.text-search[open] > summary::after { content:' \2212'; }
.facet-options { position:absolute; z-index:4; top:calc(100% + 5px); left:0; width:max-content;
          max-width:min(340px,85vw); padding:8px; border:1px solid var(--rule); border-radius:6px;
          background:var(--surface); box-shadow:0 5px 18px rgba(0,0,0,.14); display:grid; gap:3px; }
.facet-options label { display:flex; gap:7px; align-items:baseline; padding:3px 5px;
          border-radius:3px; font-size:.84rem; cursor:pointer; }
.facet-options label:hover { background:var(--surface-2); }
.facet-options small { color:var(--ink-muted); font-variant-numeric:tabular-nums; }
.more-facets { margin:0; }
.more-facets > summary, .text-search > summary { display:inline-block; }
.more-grid { display:flex; flex-wrap:wrap; gap:5px; margin-top:5px; }
.finder-foot { display:flex; flex-wrap:wrap; align-items:flex-start; gap:6px; margin-top:5px; }
.clear { font:inherit; font-size:.8rem; color:var(--accent); border:0; background:none;
         padding:2px; cursor:pointer; }
.text-search[open] { flex:1 1 300px; }
.search { display:inline; margin-left:5px; }
.search input[type=search] { font:inherit; font-size:.82rem; width:min(360px,100%); box-sizing:border-box;
          padding:3px 7px; border:1px solid var(--key-rule); border-radius:5px;
          background:var(--bg); color:var(--ink); }
.count { color:var(--ink-muted); font-size:.8rem; margin-left:auto; padding:3px 0;
         font-variant-numeric:tabular-nums; }
.protocol-table-wrap { overflow-x:auto; margin:.35em 0; border:1px solid var(--rule);
                       border-radius:5px; background:var(--surface); }
.plist { width:100%; border-collapse:collapse; font-size:.86rem; line-height:1.35; }
.plist th, .plist td { padding:7px 10px; text-align:left; vertical-align:top;
                       border-bottom:1px solid var(--rule); }
.plist thead th { position:sticky; top:0; z-index:1; background:var(--surface-2);
                  color:var(--ink-muted); font-size:.74rem; letter-spacing:.035em;
                  text-transform:uppercase; white-space:nowrap; }
.sort-button { appearance:none; border:0; padding:0; background:none; color:inherit;
               font:inherit; letter-spacing:inherit; text-transform:inherit; cursor:pointer; }
.sort-button:hover, .sort-button:focus-visible { color:var(--ink); }
.plist thead th[aria-sort="ascending"] .sort-button::after { content:' \\2191'; color:var(--accent); }
.plist thead th[aria-sort="descending"] .sort-button::after { content:' \\2193'; color:var(--accent); }
.plist tbody tr:last-child > * { border-bottom:0; }
.plist tbody tr:hover > * { background:var(--bg); }
.plist tbody th { min-width:14rem; font-weight:600; }
.plist td { color:var(--ink-muted); }
.plist .reference-cell { min-width:13rem; }
.plist .reference-cell .company { color:var(--ink-muted); }
.plist .reference-cell .no-reference { color:var(--ink-muted); }
.plist .citation-cell { text-align:right; white-space:nowrap; font-variant-numeric:tabular-nums; }
.pr .cites { color:var(--ink-muted); font-variant-numeric:tabular-nums; }
.pr .no-count { font-style:italic; }
mark { background:var(--note-bg); color:var(--ink); box-shadow:0 0 0 1px var(--note-rule);
       border-radius:2px; }
.empty { color:var(--ink-muted); font-style:italic; }
.concept { margin:1.1em 0; }
.concept p { margin:.2em 0 0; max-width:76ch; font-size:.95rem; }
.citation-filter { display:flex; flex-wrap:wrap; align-items:center; gap:5px 9px;
                   margin-top:6px; font-size:.8rem; color:var(--ink-muted); }
.citation-filter input[type=range] { width:min(240px,52vw); accent-color:var(--accent); }
.citation-filter output { min-width:7.5em; color:var(--ink); font-variant-numeric:tabular-nums; }
.citation-filter .unknown { display:flex; align-items:center; gap:5px; }
.citation-filter .snapshot { color:var(--ink-muted); }
@media (max-width:600px) {
  .card { padding:10px 12px; }
  .plist th, .plist td { padding:6px 8px; }
  .count { margin-left:0; width:100%; }
  .finder-head { display:block; }
  .facet-options { position:fixed; left:12px; right:12px; top:18%; width:auto; max-height:65vh;
                   overflow:auto; }
}
</style>"""

SEARCH_JS = r"""<script>
(function () {
  var data = JSON.parse(document.getElementById('search-data').textContent);
  var list = document.getElementById('plist');
  var items = Array.prototype.slice.call(list.children);
  var q = document.getElementById('q'), count = document.getElementById('count');
  var empty = document.getElementById('empty'), clear = document.getElementById('clear');
  var boxes = Array.prototype.slice.call(document.querySelectorAll('.facet-options input'));
  var citation = document.getElementById('citation-min');
  var citationValue = document.getElementById('citation-value');
  var keepUnknown = document.getElementById('citation-unknown');
  var citationSteps = JSON.parse(citation.dataset.steps);
  var sortButtons = Array.prototype.slice.call(document.querySelectorAll('.sort-button'));
  // Protocol name is the explicit default, rather than relying on catalogue row order.
  // This also lets showSortState() expose the initial ordering in the table header.
  var sortKey = 'name', sortDirection = 'asc';
  function norm(s) {
    return (s || '').toLowerCase().normalize('NFD').replace(/[̀-ͯ]/g, '');
  }
  function squash(s) { return s.replace(/[^a-z0-9]+/g, ''); }
  // pre-normalise every searchable field once
  var recs = data.map(function (r) {
    var f = {n: norm(r.n), a: norm(r.a.join(' | ')), c: norm(r.c + ' ' + r.m + ' ' + r.s),
             d: norm(r.d + ' ' + r.y + ' ' + r.r), t: norm(r.t)};
    f.ns = squash(f.n + ' ' + f.a);
    return f;
  });
  var W = {n: 12, a: 6, c: 3, d: 3, t: 2};
  function score(f, terms) {
    var total = 0;
    for (var i = 0; i < terms.length; i++) {
      var t = terms[i], best = 0;
      for (var k in W) if (f[k].indexOf(t) >= 0) best = Math.max(best, W[k]);
      if (!best && t.length > 2 && f.ns.indexOf(squash(t)) >= 0) best = W.a;
      if (!best) return 0;
      if (f.n.indexOf(t) === 0) best += 6;
      total += best;
    }
    return total;
  }
  function sortValue(r, key) {
    if (key === 'name') return r.n;
    if (key === 'reference') return r.r;
    if (key === 'citations') return r.x;
    if (key === 'family') return r.p.assay.join(', ');
    if (key === 'scale') return r.p.scale.join(', ');
    if (key === 'platform') return r.p.platform.join(', ');
    return '';
  }
  function orderHits(hits) {
    if (!sortKey) {
      hits.sort(function (x, y) { return y[0] - x[0] || x[1] - y[1]; });
      return;
    }
    var sign = sortDirection === 'asc' ? 1 : -1;
    hits.sort(function (x, y) {
      var a = sortValue(data[x[1]], sortKey), b = sortValue(data[y[1]], sortKey);
      // Unknown citation counts are metadata gaps, not zero; keep them after known
      // values in either direction.
      if (sortKey === 'citations') {
        if (a === null || b === null) return a === b ? x[1] - y[1] : (a === null ? 1 : -1);
        if (a !== b) return sign * (a - b);
      } else {
        var compared = String(a).localeCompare(String(b), undefined,
                                                {numeric:true, sensitivity:'base'});
        if (compared) return sign * compared;
      }
      return data[x[1]].n.localeCompare(data[y[1]].n, undefined,
                                        {numeric:true, sensitivity:'base'});
    });
  }
  function showSortState() {
    sortButtons.forEach(function (button) {
      var active = button.dataset.sort === sortKey;
      button.parentNode.setAttribute('aria-sort', active
        ? (sortDirection === 'asc' ? 'ascending' : 'descending') : 'none');
    });
  }
  function run() {
    var terms = norm(q.value).split(/\s+/).filter(Boolean);
    var active = {};
    boxes.forEach(function (box) {
      if (box.checked) (active[box.dataset.facet] || (active[box.dataset.facet] = [])).push(box.value);
    });
    var hits = [];
    var minCitations = citationSteps[Number(citation.value)];
    citationValue.textContent = minCitations ? minCitations.toLocaleString() + '+' : 'any';
    for (var i = 0; i < recs.length; i++) {
      var r = data[i], f = recs[i];
      var matches = Object.keys(active).every(function (facet) {
        return active[facet].some(function (value) { return r.p[facet].indexOf(value) >= 0; });
      });
      if (!matches) continue;
      if (r.x === null ? !keepUnknown.checked : r.x < minCitations) continue;
      var s = terms.length ? score(f, terms) : 1;
      if (s) hits.push([s, i]);
    }
    orderHits(hits);
    var shown = {};
    hits.forEach(function (h) {
      var li = items[h[1]];
      shown[h[1]] = 1;
      list.appendChild(li);
      li.hidden = false;
    });
    items.forEach(function (li, i) { if (!shown[i]) li.hidden = true; });
    count.textContent = hits.length + ' of ' + items.length + ' protocols';
    empty.hidden = hits.length > 0;
    clear.hidden = !q.value && citation.value === '0' && keepUnknown.checked
      && !boxes.some(function (box) { return box.checked; });
    document.querySelectorAll('.facet').forEach(function (group) {
      var n = group.querySelectorAll('input:checked').length, summary = group.querySelector('summary');
      summary.classList.toggle('active', n > 0);
      summary.querySelector('.selected').textContent = n ? ' (' + n + ')' : '';
    });
    var h = [];
    if (q.value) h.push('q=' + encodeURIComponent(q.value));
    if (minCitations) h.push('c=' + minCitations);
    if (!keepUnknown.checked) h.push('unknown=0');
    if (sortKey) {
      h.push('sort=' + encodeURIComponent(sortKey));
      h.push('order=' + sortDirection);
    }
    boxes.forEach(function (box) {
      if (box.checked) h.push('f=' + encodeURIComponent(box.dataset.facet + '~' + box.value));
    });
    try { history.replaceState(null, '', h.length ? '#' + h.join('&') : location.pathname); }
    catch (err) { /* file:// in some browsers */ }
  }
  q.addEventListener('input', run);
  citation.addEventListener('input', run);
  keepUnknown.addEventListener('change', run);
  boxes.forEach(function (box) { box.addEventListener('change', run); });
  sortButtons.forEach(function (button) {
    button.addEventListener('click', function () {
      if (sortKey === button.dataset.sort) {
        sortDirection = sortDirection === 'asc' ? 'desc' : 'asc';
      } else {
        sortKey = button.dataset.sort;
        sortDirection = button.dataset.first || 'asc';
      }
      showSortState();
      run();
    });
  });
  var groups = Array.prototype.slice.call(document.querySelectorAll('.facet'));
  groups.forEach(function (group) {
    group.addEventListener('toggle', function () {
      if (group.open) groups.forEach(function (other) { if (other !== group) other.open = false; });
    });
  });
  document.addEventListener('click', function (ev) {
    if (!ev.target.closest('.facet')) groups.forEach(function (group) { group.open = false; });
  });
  clear.addEventListener('click', function () {
    q.value = ''; citation.value = '0'; keepUnknown.checked = true;
    boxes.forEach(function (box) { box.checked = false; }); run();
  });
  q.addEventListener('keydown', function (ev) { if (ev.key === 'Escape') { q.value = ''; run(); } });
  document.addEventListener('keydown', function (ev) {
    if (ev.key === '/' && document.activeElement !== q) { ev.preventDefault(); q.focus(); }
    if (ev.key === 'Escape') groups.forEach(function (group) { group.open = false; });
  });
  location.hash.replace(/^#/, '').split('&').forEach(function (kv) {
    var p = kv.split('='), v = decodeURIComponent(p[1] || '');
    if (p[0] === 'q') q.value = v;
    if (p[0] === 'c') {
      var wanted = Number(v), at = citationSteps.indexOf(wanted);
      if (at >= 0) citation.value = String(at);
    }
    if (p[0] === 'unknown' && v === '0') keepUnknown.checked = false;
    if (p[0] === 'sort' && sortButtons.some(function (button) {
      return button.dataset.sort === v;
    })) sortKey = v;
    if (p[0] === 'order' && (v === 'asc' || v === 'desc')) sortDirection = v;
    if (p[0] === 'f') {
      var bits = v.split('~'), facet = bits.shift(), value = bits.join('~');
      boxes.forEach(function (box) {
        if (box.dataset.facet === facet && box.value === value) box.checked = true;
      });
    }
  });
  showSortState();
  run();
})();
</script>"""


PRIMARY_FACETS = ("scale", "assay", "platform", "index_introduction", "index_architecture",
                  "partitioning", "amplification")


def facet_group(key: str, counts: dict[str, dict[str, int]]) -> str:
    label, allowed = props.FILTER_FACETS[key]
    options = []
    for i, value in enumerate(allowed):
        n = counts[key][value]
        if not n:
            continue
        ident = f"f-{key}-{i}"
        options.append(
            f'<label for="{ident}"><input id="{ident}" type="checkbox" '
            f'data-facet="{e(key)}" value="{e(value)}">'
            f'<span>{e(value)}</span> <small>{n}</small></label>')
    return (f'<details class="facet"><summary>{e(label)}'
            f'<span class="selected"></span></summary>'
            f'<div class="facet-options">{"".join(options)}</div></details>')


def facet_controls(ps: list[dict]) -> str:
    counts = props.option_counts(ps)
    primary = "".join(facet_group(key, counts) for key in PRIMARY_FACETS)
    more = "".join(facet_group(key, counts) for key in props.FILTER_FACETS
                   if key not in PRIMARY_FACETS)
    return f"""<div class="facetbar">
{primary}
<details class="more-facets"><summary>More filters</summary>
<div class="more-grid">{more}</div></details>
</div>"""


def citation_controls(ps: list[dict]) -> str:
    """A useful nonlinear citation threshold plus explicit handling of absent data."""
    maximum = max((p["citations"] for p in ps if p["citations"] is not None), default=0)
    steps = [x for x in (0, 10, 25, 50, 100, 250, 500, 1000, 2500, 5000, 10000)
             if x <= maximum]
    if maximum and maximum not in steps:
        steps.append(maximum)
    encoded = e(json.dumps(steps, separators=(",", ":")), quote=True)
    dates = sorted({p["citations_retrieved"] for p in ps if p["citations_retrieved"]})
    dated = dates[0] if len(dates) == 1 else "various dates"
    return f"""<div class="citation-filter">
<label for="citation-min">Minimum citations</label>
<input id="citation-min" type="range" min="0" max="{len(steps) - 1}" value="0"
       step="1" data-steps="{encoded}">
<output id="citation-value" for="citation-min">any</output>
<label class="unknown"><input id="citation-unknown" type="checkbox" checked>
Keep protocols without citation data</label>
<span class="snapshot">OpenAlex snapshot · {e(dated)}</span>
</div>"""


def search_section(ps: list[dict]) -> str:
    items = "\n".join(result_item(i, p) for i, p in enumerate(ps))
    return f"""<form class="finder" role="search" onsubmit="return false">
<div class="finder-head"><h2>Filter protocols</h2></div>
{facet_controls(ps)}
{citation_controls(ps)}
<div class="finder-foot"><details class="text-search"><summary>Text search</summary>
<span class="search"><label for="q" class="sr-only">Search protocol text</label>
<input id="q" type="search" placeholder="Name, paper or description"
       autocomplete="off" spellcheck="false"></span></details>
<span class="count" id="count" aria-live="polite">{len(ps)} protocols</span>
<button class="clear" id="clear" type="button" hidden>Clear filters</button></div>
</form>
<div class="protocol-table-wrap"><table class="plist">
<thead><tr><th scope="col" aria-sort="none"><button type="button" class="sort-button"
data-sort="name">Protocol</button></th>
<th scope="col" aria-sort="none"><button type="button" class="sort-button"
data-sort="reference">Paper / company</button></th>
<th scope="col" aria-sort="none"><button type="button" class="sort-button"
data-sort="citations" data-first="desc">Citations</button></th>
<th scope="col" aria-sort="none"><button type="button" class="sort-button"
data-sort="family">Protocol family</button></th>
<th scope="col" aria-sort="none"><button type="button" class="sort-button"
data-sort="scale">Bulk / single-cell</button></th>
<th scope="col" aria-sort="none"><button type="button" class="sort-button"
data-sort="platform">Target sequencer</button></th></tr></thead>
<tbody id="plist">{items}</tbody>
</table></div>
<p class="empty" id="empty" hidden>Nothing matches. Try fewer words, or clear the filters.</p>
<script type="application/json" id="search-data">{search_data(ps)}</script>
{SEARCH_JS}"""


def build(omit: set[str] = frozenset(), run_checks: bool = True) -> tuple[str, int]:
    every = collect(omit, run_checks)
    # Keep the complete published catalogue visible. Reference notes and work in progress
    # remain out of the reader-facing index; entries without a finished schematic are
    # listed by name but do not send readers into the notes.
    # Keep the non-JavaScript table in the same explicit A-Z order as the interactive
    # default. JavaScript may subsequently apply a sort encoded in the URL.
    pub = sorted((p for p in every if p["section"] == "published"),
                 key=lambda p: p["name"].casefold())
    total = (sum(p["checks"] for p in every)
             + (sum(checks(d)[0] for d in EXTRA_SUITES) if run_checks else 0))
    body = f"""<div class="wrap">
{search_section(pub)}
</div>"""
    page = ("<!doctype html>\n<html lang=\"en\">\n<head>\n"
            + head("NGS 101", include_header=False)
            + "\n" + EXTRA_CSS + "\n</head>\n<body class=\"index\">\n"
            + site_header("index.html") + "\n" + body
            + "\n</body>\n</html>\n")
    return page, total


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--omit", action="append", default=[], metavar="DIR",
                    help="list DIR without its diagram page (its build failed)")
    ap.add_argument("--no-checks", action="store_true",
                    help="do not run the self-tests for the check counts")
    ap.add_argument("-o", "--out", type=Path, default=OUT)
    args = ap.parse_args(argv)
    page, total = build(set(args.omit), not args.no_checks)
    args.out.write_text(page, encoding="utf-8")
    print(f"wrote {args.out}  ({args.out.stat().st_size:,} bytes)  {total} checks total")
    return 0


if __name__ == "__main__":
    sys.exit(main())
