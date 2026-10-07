#!/usr/bin/env python3
"""Generate index.html -- the landing page, with a searchable list of the protocols.

Nothing about an individual protocol is typed here. The list comes from
`catalogue/ours.tsv` (one row per directory; its `section` column says whether it is a
published protocol or our own work in progress), joined to the scraped catalogue for
category, family members, papers and year, and to the notes themselves: a protocol's
title is its first note's `# ` heading and its blurb that note's first real paragraph.
So a new `<slug>/01_*.md` plus a row in ours.tsv is all it takes to appear here.

The published protocols get a client-side search: every field above plus the full text of
the notes goes into a JSON index embedded in the page, and a few lines of vanilla JS
filter the pre-rendered list. No library, no fetch() -- it works from file:// too, and
without JavaScript the full list is simply shown.

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
from mdfacts import expand  # noqa: E402
from mdrender import render  # noqa: E402
from page import caveat, head, info, legend  # noqa: E402

OUT = ROOT / "index.html"

# Suites that are not a protocol directory but still count towards the total.
EXTRA_SUITES = ["catalogue", "gcbias", "gcbias/datasets"]

BLURB_CHARS = 330          # a blurb is cut at a sentence end before this many characters
OUT_RE = re.compile(r'^OUT\s*=\s*HERE\.parent\s*/\s*"([^"]+\.html)"', re.M)


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
    """The page a protocol's tools/build_page.py writes (its `OUT = HERE.parent / ...`)."""
    bp = ROOT / d / "tools" / "build_page.py"
    if not bp.exists():
        return None
    m = OUT_RE.search(bp.read_text(encoding="utf-8"))
    return ROOT / d / m.group(1) if m else None


def checks(d: str) -> tuple[int, bool]:
    """Run a directory's self-test and report (count, passed); (0, False) if none."""
    st = ROOT / d / "tools" / "selftest.py"
    if not st.exists():
        return 0, False
    r = subprocess.run([sys.executable, str(st)], capture_output=True, text=True)
    m = re.search(r"all (\d+) checks passed", r.stdout)
    if m:
        return int(m.group(1)), True
    m = re.search(r"(\d+) of (\d+) checks FAILED", r.stdout)
    return (int(m.group(2)) if m else 0), False


# ------------------------------------------------------------- one protocol
CATEGORY_LABEL = {"TODO list": "on the scg_lib_structs TODO list",
                  "Ours, not in scg_lib_structs": "not in scg_lib_structs"}


def entry(r: dict, omit: set[str], run_checks: bool) -> dict:
    d = r["dir"]
    crows = cat.rows_for_dir(d)
    defining = next((x for x in crows if x["is_defining"] == "yes"), {})
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
    has_page = bool(page and page.exists() and d not in omit)
    n, ok = checks(d) if run_checks else (0, False)
    category = crows[0]["category"] if crows else ""
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
    }


def collect(omit: set[str] = frozenset(), run_checks: bool = True) -> list[dict]:
    return [entry(r, set(omit), run_checks) for r in cat.ours()]


# ------------------------------------------------------------------- render
e = html.escape


def status_label(p: dict) -> str:
    return "notes + diagram page" if p["page"] else "notes"


def links(p: dict) -> str:
    out = []
    if p["page"]:
        out.append(f'<a class="go" href="{e(p["page"])}">diagram page &rarr;</a>')
    for n in p["notes"]:
        out.append(f'<a href="{e(n["href"])}">{e(n["title"])}</a>')
    if p["scg_page"]:
        out.append(f'<a class="ext" href="{e(p["scg_page"])}">scg_lib_structs</a>')
    return '<div class="links">' + "".join(out) + "</div>"


def meta_line(p: dict) -> str:
    bits = []
    if p["year"]:
        bits.append(e(p["year"]))
    if p["doi"]:
        bits.append(f'<a href="https://doi.org/{e(p["doi"])}">doi:{e(p["doi"])}</a>')
    bits.append(f"<code>{e(p['dir'])}/</code>")
    if p["checks"]:
        bits.append(f"{p['checks']} checks passing" if p["checks_ok"]
                    else f"&#9888; {p['checks']} checks, suite FAILING")
    return '<div class="meta">' + " &middot; ".join(bits) + "</div>"


def badges(p: dict) -> str:
    b = [f'<span class="tag mod">{e(p["modality"])}</span>']
    if p["category"]:
        b.append(f'<span class="tag">{e(p["category"])}</span>')
    b.append(f'<span class="tag{" pg" if p["page"] else ""}">{status_label(p)}</span>')
    return '<div class="tags">' + "".join(b) + "</div>"


def result_item(i: int, p: dict) -> str:
    first = p["notes"][0]["href"] if p["notes"] else (p["page"] or "")
    name = f'<a href="{e(first)}">{e(p["name"])}</a>' if first else e(p["name"])
    aka = (f'<p class="aka">also: {e(", ".join(p["aka"]))}</p>' if p["aka"] else "")
    return f"""<li class="pr" data-i="{i}" id="p-{e(p['dir'])}">
<h3>{name}</h3>
{badges(p)}
{aka}<p class="blurb">{e(p['blurb'])}</p>
<p class="snip" hidden></p>
{links(p)}
{meta_line(p)}
</li>"""


def wip_card(p: dict) -> str:
    return f"""<div class="card">
<span class="kind">{e(p['modality'])} &middot; work in progress &middot; status {e(p['status'])}</span>
<h3>{e(p['name'])}</h3>
<p>{e(p['blurb'])}</p>
{links(p)}
{meta_line(p)}
</div>"""


def search_data(ps: list[dict]) -> str:
    """The search index: one record per published protocol, in list order."""
    recs = [{"n": p["name"], "a": p["aka"], "c": p["category"], "m": p["modality"],
             "s": f'{p["status"]} {status_label(p)}', "pg": bool(p["page"]),
             "d": p["doi"], "y": p["year"],
             "t": " ".join([*p["papers"], *(n["title"] for n in p["notes"])]),
             "x": p["text"]} for p in ps]
    # embedded in a <script>, so no "</" may appear literally
    return json.dumps(recs, ensure_ascii=False, separators=(",", ":")).replace("</", "<\\/")


def concepts() -> str:
    out = []
    for f in sorted((ROOT / "ref" / "concepts").glob("*.md")):
        n = read_note(f)
        out.append(f'<div class="concept"><h3><a href="{e(n["href"])}">{e(n["title"])}</a>'
                   f'</h3><p>{e(n["blurb"])}</p></div>')
    return "\n".join(out)


EXTRA_CSS = """<style>
.cards { display:grid; gap:14px; grid-template-columns:repeat(auto-fit,minmax(min(330px,100%),1fr));
         margin:1.2em 0; align-items:start; }
.card { border:1px solid var(--rule); border-radius:6px; padding:16px 18px;
        background:var(--surface); display:flex; flex-direction:column; gap:8px; }
.card h3, .pr h3, .concept h3 { margin:0; font-size:1.06rem; }
.card h3::before, .pr h3::before, .concept h3::before { content:none; }
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

/* search */
.search { position:sticky; top:0; z-index:2; background:var(--bg); padding:10px 0 8px;
          display:flex; flex-direction:column; gap:8px; border-bottom:1px solid var(--rule); }
.search input[type=search] { font:inherit; font-size:1rem; width:100%; box-sizing:border-box;
          padding:9px 12px; border:1px solid var(--key-rule); border-radius:6px;
          background:var(--surface); color:var(--ink); }
.facets { display:flex; flex-wrap:wrap; gap:6px; align-items:center; font-size:.85rem; }
.facets button { font:inherit; font-size:.82rem; padding:3px 10px; border-radius:999px;
          border:1px solid var(--rule); background:var(--surface); color:var(--ink-muted);
          cursor:pointer; }
.facets button[aria-pressed=true] { background:var(--accent); border-color:var(--accent);
          color:var(--bg); }
.facets .sep { width:1px; height:1.2em; background:var(--rule); margin:0 4px; }
.count { color:var(--ink-muted); font-size:.85rem; margin-left:auto;
         font-variant-numeric:tabular-nums; }
.plist { list-style:none; padding:0; margin:.6em 0; display:grid; gap:10px; }
.pr { border:1px solid var(--rule); border-radius:6px; padding:12px 16px;
      background:var(--surface); display:flex; flex-direction:column; gap:6px; }
.pr p { margin:0; font-size:.93rem; max-width:90ch; }
.pr .aka { font-size:.84rem; color:var(--ink-muted); }
.pr .snip { font-size:.84rem; color:var(--ink-muted); border-left:2px solid var(--key-rule);
            padding-left:8px; }
mark { background:var(--note-bg); color:var(--ink); box-shadow:0 0 0 1px var(--note-rule);
       border-radius:2px; }
.empty { color:var(--ink-muted); font-style:italic; }
.concept { margin:1.1em 0; }
.concept p { margin:.2em 0 0; max-width:76ch; font-size:.95rem; }
@media (max-width:600px) {
  .search { position:static; }      /* a sticky search box would eat a phone screen */
  .pr, .card { padding:10px 12px; }
  .count { margin-left:0; width:100%; }
}
</style>"""

SEARCH_JS = r"""<script>
(function () {
  var data = JSON.parse(document.getElementById('search-data').textContent);
  var list = document.getElementById('plist');
  var items = Array.prototype.slice.call(list.children);
  var q = document.getElementById('q'), count = document.getElementById('count');
  var empty = document.getElementById('empty');
  var facets = document.querySelectorAll('.facets button');
  var state = {m: '', pg: ''};
  function norm(s) {
    return (s || '').toLowerCase().normalize('NFD').replace(/[̀-ͯ]/g, '');
  }
  function squash(s) { return s.replace(/[^a-z0-9]+/g, ''); }
  // pre-normalise every searchable field once
  var recs = data.map(function (r) {
    var f = {n: norm(r.n), a: norm(r.a.join(' | ')), c: norm(r.c + ' ' + r.m + ' ' + r.s),
             d: norm(r.d + ' ' + r.y), t: norm(r.t), x: norm(r.x)};
    f.ns = squash(f.n + ' ' + f.a);
    return f;
  });
  var W = {n: 12, a: 6, c: 3, d: 3, t: 2, x: 1};
  function score(f, terms) {
    var total = 0;
    for (var i = 0; i < terms.length; i++) {
      var t = terms[i], best = 0;
      for (var k in W) if (f[k].indexOf(t) >= 0) best = Math.max(best, W[k]);
      if (!best && t.length > 2 && f.ns.indexOf(squash(t)) >= 0) best = W.a;
      if (!best) return 0;
      if (best === W.x) best += Math.min(f.x.split(t).length - 1, 10) / 10;
      if (f.n.indexOf(t) === 0) best += 6;
      total += best;
    }
    return total;
  }
  function snippet(raw, f, terms, phrase) {
    // only when the match is in the note text, not already visible in the card;
    // the whole query as a phrase if it occurs, else the earliest single term
    var pos = -1, len = 0;
    if (terms.length > 1 && f.n.indexOf(phrase) < 0) {
      pos = f.x.indexOf(phrase); len = phrase.length;
    }
    if (pos < 0) for (var i = 0; i < terms.length; i++) {
      var t = terms[i];
      if (f.n.indexOf(t) >= 0 || f.a.indexOf(t) >= 0 || f.c.indexOf(t) >= 0 ||
          f.d.indexOf(t) >= 0) continue;
      var p = f.x.indexOf(t);
      if (p >= 0 && (pos < 0 || p < pos)) { pos = p; len = t.length; }
    }
    if (pos < 0) return null;
    var a = Math.max(0, pos - 90), b = Math.min(raw.length, pos + len + 110);
    return [(a > 0 ? '…' : '') + raw.slice(a, pos), raw.slice(pos, pos + len),
            raw.slice(pos + len, b) + (b < raw.length ? '…' : '')];
  }
  function run() {
    var terms = norm(q.value).split(/\s+/).filter(Boolean);
    var hits = [];
    for (var i = 0; i < recs.length; i++) {
      var r = data[i], f = recs[i];
      if (state.m && r.m !== state.m) continue;
      if (state.pg === 'page' && !r.pg) continue;
      if (state.pg === 'notes' && r.pg) continue;
      var s = terms.length ? score(f, terms) : 1;
      if (s) hits.push([s, i]);
    }
    hits.sort(function (x, y) { return y[0] - x[0] || x[1] - y[1]; });
    var shown = {};
    hits.forEach(function (h) {
      var li = items[h[1]], sn = li.querySelector('.snip');
      shown[h[1]] = 1;
      list.appendChild(li);
      li.hidden = false;
      var parts = terms.length && snippet(data[h[1]].x, recs[h[1]], terms, terms.join(' '));
      // norm() keeps offsets (NFD, then the combining marks are dropped again), so the
      // match position in the normalised text is the position in the raw text
      if (parts) {
        sn.textContent = '';
        sn.appendChild(document.createTextNode(parts[0]));
        var m = document.createElement('mark'); m.textContent = parts[1];
        sn.appendChild(m);
        sn.appendChild(document.createTextNode(parts[2]));
        sn.hidden = false;
      } else { sn.hidden = true; }
    });
    items.forEach(function (li, i) { if (!shown[i]) li.hidden = true; });
    count.textContent = hits.length + ' of ' + items.length + ' protocols';
    empty.hidden = hits.length > 0;
    var h = [];
    if (q.value) h.push('q=' + encodeURIComponent(q.value));
    if (state.m) h.push('m=' + state.m);
    if (state.pg) h.push('pg=' + state.pg);
    try { history.replaceState(null, '', h.length ? '#' + h.join('&') : location.pathname); }
    catch (err) { /* file:// in some browsers */ }
  }
  function press() {
    facets.forEach(function (b) {
      b.setAttribute('aria-pressed', String(state[b.dataset.k] === b.dataset.v));
    });
  }
  facets.forEach(function (b) {
    b.addEventListener('click', function () { state[b.dataset.k] = b.dataset.v; press(); run(); });
  });
  q.addEventListener('input', run);
  q.addEventListener('keydown', function (ev) { if (ev.key === 'Escape') { q.value = ''; run(); } });
  document.addEventListener('keydown', function (ev) {
    if (ev.key === '/' && document.activeElement !== q) { ev.preventDefault(); q.focus(); }
  });
  location.hash.replace(/^#/, '').split('&').forEach(function (kv) {
    var p = kv.split('='), v = decodeURIComponent(p[1] || '');
    if (p[0] === 'q') q.value = v; else if (p[0] === 'm' || p[0] === 'pg') state[p[0]] = v;
  });
  document.getElementById('facets').hidden = false;
  press(); run();
})();
</script>"""


def search_section(ps: list[dict]) -> str:
    items = "\n".join(result_item(i, p) for i, p in enumerate(ps))
    mods = [m for m in cat.MODALITIES if any(p["modality"] == m for p in ps)]
    mbtn = "".join(f'<button type="button" data-k="m" data-v="{m}">{m}</button>' for m in mods)
    return f"""<form class="search" role="search" onsubmit="return false">
<label for="q" class="kind">Search {len(ps)} published protocols &mdash; name, family member,
category, DNA/RNA, DOI, year, or anything in the notes (press <kbd>/</kbd>)</label>
<input id="q" type="search" placeholder="e.g. template switching, Tn5, 10x, UMI, 2017, nbt.2282"
       autocomplete="off" spellcheck="false">
<div class="facets" id="facets" hidden>
<button type="button" data-k="m" data-v="">all</button>{mbtn}
<span class="sep"></span>
<button type="button" data-k="pg" data-v="">any status</button><button type="button" data-k="pg" data-v="page">with diagram page</button><button type="button" data-k="pg" data-v="notes">notes only</button>
<span class="count" id="count" aria-live="polite">{len(ps)} protocols</span>
</div>
</form>
<ol class="plist" id="plist">
{items}
</ol>
<p class="empty" id="empty" hidden>Nothing matches. Try fewer words, or clear the filters.</p>
<script type="application/json" id="search-data">{search_data(ps)}</script>
{SEARCH_JS}"""


def build(omit: set[str] = frozenset(), run_checks: bool = True) -> tuple[str, int]:
    every = collect(omit, run_checks)
    pub = [p for p in every if p["section"] == "published"]
    wip = [p for p in every if p["section"] == "wip"]
    # protocols with a diagram page first, then the rest in ours.tsv order
    pub.sort(key=lambda p: not p["page"])
    total = (sum(p["checks"] for p in every)
             + (sum(checks(d)[0] for d in EXTRA_SUITES) if run_checks else 0))
    n_notes = len([q for q in ROOT.rglob("*.md")
                   if not {".git", "_data", "pdf", "__pycache__", ".cache", "_site"}
                   & set(q.relative_to(ROOT).parts)])
    n_pages = sum(1 for p in pub if p["page"])
    checks_line = (f" &mdash; <b>{total} checks</b> across the protocols" if total else "")
    body = f"""<div class="wrap">
<h1>NGS protocol chemistry</h1>

{info("""What the DNA actually looks like at every step of a sequencing protocol, drawn
base-by-base in the idiom of the Teichmann lab's
<a href="https://teichlab.github.io/scg_lib_structs/">scg_lib_structs</a> pages.""")}

{legend(f"""<b>Everything here is generated, not transcribed.</b> Each protocol's construct
is defined once as a table of segments; the diagrams, oligo lists and final structures are
derived from it, and a self-test asserts the facts that were checked by hand{checks_line}.
Where a real plasmid exists, primer sites and amplicon sizes are computed against the
actual map rather than quoted from a paper.""")}

<h2 id="protocols">Published protocols</h2>
{info(f"""{len(pub)} published methods, each with reference notes read from its own
paper(s) &mdash; oligos verbatim, how they interlock, step by step, and the final library.
{n_pages} so far also have a generated base-by-base diagram page.""")}
{search_section(pub)}

<h2 id="wip">Work in progress</h2>
<div class="wipbox">
{caveat("""<b>Our own, unpublished or still-moving designs.</b> These are not published
protocols: the chemistry is being worked out here, the notes change, and some of the final
library is predicted rather than documented. Treat them as lab notebooks.""")}
<div class="cards">
{"".join(wip_card(p) for p in wip)}
</div>
</div>

<h2 id="catalogue">The catalogue &mdash; what is out there</h2>
{info(f"""Documenting a protocol starts with knowing it exists. <b>{cat.n_protocols()}
protocols</b> and <b>{cat.n_papers()} papers</b> are catalogued in
<code>catalogue/scg_lib_structs.tsv</code>, scraped from
<a href="https://github.com/Teichlab/scg_lib_structs">Teichlab/scg_lib_structs</a> and
resolved against Crossref and NCBI &mdash; {cat.n_documented()} with a drawn page upstream
to check ourselves against, {cat.n_todo()} named but never drawn. {cat.n_ours()} of them
have a directory here. One row per (protocol, paper), because a method can have several
papers and a paper can define several methods.
<a href="catalogue/README.html"><b>read the worklist &rarr;</b></a>""")}

<h2 id="concepts">Chemistry that recurs</h2>
{info("""Most new protocols are a new front end bolted onto an old back end. These
pieces turn up again and again, so they are written once and shared rather than repeated
per protocol.""")}
{concepts()}

<h2 id="notes">Written notes</h2>
{info(f"""Alongside the diagram pages, the repository carries {n_notes} Markdown notes
&mdash; protocol readings, dataset inventories, the enzyme/buffer comparison, and working
logs. All of them are rendered to HTML by <code>build_docs.py</code>:
<a href="notes.html"><b>browse all notes &rarr;</b></a>""")}

<h2>How it is built</h2>
{info("""<code>lib/</code> holds the protocol-agnostic machinery: a renderer for
character-aligned duplex diagrams, the canonical Illumina and NEBNext sequences, Tn5 and
reverse-transcription building blocks, a GenBank parser with restriction digestion and
circular-aware PCR prediction, and the padlock capture model. Each protocol directory holds
its own notes, its segment definitions and its self-test, which runs the shared checks plus
its own. The Markdown notes go through <code>lib/mdrender.py</code>, a small
dependency-free renderer, driven by <code>build_docs.py</code>; this page and the search
index come from <code>build_index.py</code>, and <code>build_site.py</code> assembles the
lot into the published site.""")}

{caveat("""<b>On trusting this.</b> The checks are regression tests, not proof. The useful
ones assert relationships that <i>must</i> hold if a model is right &mdash; for instance the
Atrandi read layout is derived from the segment table and then asserted against offsets
hard-coded in a demultiplexer that was written from sequencing data. When a model reproduces
a number nobody fed it, that is worth locking down. Anything inferred rather than documented
is marked as such on the page it appears on.""")}
</div>"""
    page = ("<!doctype html>\n<html lang=\"en\">\n<head>\n" + head("NGS Protocol Chemistry")
            + "\n" + EXTRA_CSS + "\n</head>\n<body>\n" + body + "\n</body>\n</html>\n")
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
