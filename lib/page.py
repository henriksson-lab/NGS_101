"""
Shared page scaffolding for protocol chemistry pages.

The stylesheet lives here rather than in each protocol's build script so that every page
looks like the same publication and a fix reaches all of them. Region colours are tokens,
defined for all three theme states (bare :root, system-dark, explicitly-stamped dark);
`inf` marks how well a region is KNOWN, orthogonally to what it IS.
"""

from __future__ import annotations

import html

STYLE = """<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;600&family=IBM+Plex+Sans:wght@400;500;600&display=swap">
<style>
/* ---- tokens: complete light palette on bare :root ------------------------ */
:root {
  --bg:#fbfcfd; --surface:#f2f6fa; --surface-2:#e9eff5;
  --ink:#16202b; --ink-muted:#566879; --rule:#dce4ed;
  --link:#0b5fa5; --accent:#0b5fa5;
  --note-bg:#fff6e3; --note-rule:#d79a16;
  --key-bg:#eef3f8; --key-rule:#7f95aa;
  --c-p5:#08519c; --c-p7:#a50f15; --c-s5:#3d87bd; --c-s7:#e2663f;
  --c-me:#7d7d7d; --c-t7:#1f4fd8; --c-cbc:#d6417c; --c-umi:#6f58b5;
  --c-r1:#4a1486; --c-r2:#6a51a3; --c-r3:#807dba; --c-tso:#1f8a4d; --c-w1:#d6341a;
  --inf:#8a98a6;
  --mono:'IBM Plex Mono','DejaVu Sans Mono',Menlo,Consolas,monospace;
  --sans:'IBM Plex Sans','Helvetica Neue',Helvetica,Arial,sans-serif;
}
/* ---- same tokens, dark values, for system-dark (un-stamped) -------------- */
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    --bg:#0e1419; --surface:#161e26; --surface-2:#1d272f;
    --ink:#e2e9ef; --ink-muted:#93a4b4; --rule:#28333d;
    --link:#78b4ee; --accent:#78b4ee;
    --note-bg:#2a2214; --note-rule:#b98a1f;
    --key-bg:#17202a; --key-rule:#5a7389;
    --c-p5:#6ba6e2; --c-p7:#f0837f; --c-s5:#79bde4; --c-s7:#f59b78;
    --c-me:#9ba3ab; --c-t7:#8aa5ff; --c-cbc:#ff8ab8; --c-umi:#b0a2e8;
    --c-r1:#bda4ec; --c-r2:#b9a6e6; --c-r3:#a99ee0; --c-tso:#4fc98a; --c-w1:#ff8360;
    --inf:#78868f;
  }
}
/* ---- and again for an explicit dark choice ------------------------------ */
:root[data-theme="dark"] {
  --bg:#0e1419; --surface:#161e26; --surface-2:#1d272f;
  --ink:#e2e9ef; --ink-muted:#93a4b4; --rule:#28333d;
  --link:#78b4ee; --accent:#78b4ee;
  --note-bg:#2a2214; --note-rule:#b98a1f;
  --key-bg:#17202a; --key-rule:#5a7389;
  --c-p5:#6ba6e2; --c-p7:#f0837f; --c-s5:#79bde4; --c-s7:#f59b78;
  --c-me:#9ba3ab; --c-t7:#8aa5ff; --c-cbc:#ff8ab8; --c-umi:#b0a2e8;
  --c-r1:#bda4ec; --c-r2:#b9a6e6; --c-r3:#a99ee0; --c-tso:#4fc98a; --c-w1:#ff8360;
  --inf:#78868f;
}

body { background:var(--bg); color:var(--ink); font-family:var(--sans);
       margin:0; padding-block:32px 96px; padding-inline:16px; line-height:1.5; }
.wrap { max-width:1180px; margin:0 auto; display:flex; flex-direction:column; gap:4px; }

h1 { font-size:clamp(1.5rem,1.1rem + 1.6vw,2.1rem); font-weight:600; line-height:1.2;
     text-wrap:balance; margin:0 0 .2em; letter-spacing:-0.01em; }
.research-notes { margin:0 0 .7em; font-size:.86rem; }
.research-notes a { color:var(--ink-muted); }
h2 { font-size:1.3rem; font-weight:600; text-wrap:balance; margin:2.2em 0 .2em;
     padding-bottom:.3em; border-bottom:2px solid var(--rule); }
h3 { font-size:1rem; font-weight:600; text-wrap:balance; margin:2em 0 .3em;
     color:var(--ink); }
h3::before { content:''; display:block; width:28px; height:2px; background:var(--accent);
             margin-bottom:.6em; opacity:.55; }
p { margin:.5em 0; }
info { display:block; max-width:68ch; color:var(--ink); }
a { color:var(--link); text-underline-offset:2px; }
a:focus-visible, :focus-visible { outline:2px solid var(--accent); outline-offset:2px; }
code { font-family:var(--mono); font-size:.9em; background:var(--surface-2);
       padding:.08em .3em; border-radius:3px; }

seq { font-family:var(--mono); font-size:.86rem; display:block; margin:.4em 0 0; }
seq p { margin:.3em 0; text-indent:-2.9em; padding-left:2.9em; }

pre { overflow-x:auto; margin:.6em 0 .2em; background:var(--surface);
      border-left:3px solid var(--rule); padding:12px 14px; }
align { font-family:var(--mono); display:block; white-space:pre; }
align.long  { font-size:.74rem; line-height:1.25; }
align.small { font-size:.84rem; line-height:1.3; }
align i { color:var(--ink-muted); display:block; font-style:italic; line-height:1.5;
          margin-bottom:.7em; max-width:96ch; white-space:pre-wrap; }

/* Reaction/construct diagrams stay aligned with prose until their intrinsic content is
   wider; only then do they grow toward the viewport. Long content still scrolls. */
.chem-panel { margin:.65em 0 .25em; padding:0; min-width:0; width:100%; max-width:100%; }
.chem-panel figcaption { color:var(--ink-muted); font-style:italic; font-size:.86rem;
                         line-height:1.45; margin:0 0 .35em; max-width:96ch; }
.diagram-scroll { overflow-x:auto; overflow-y:hidden; width:max-content; min-width:100%;
                  max-width:calc(100vw - 32px); box-sizing:border-box;
                  position:relative; left:50%; transform:translateX(-50%);
                  scrollbar-gutter:stable; overscroll-behavior-inline:contain;
                  background:var(--surface); border-left:3px solid var(--rule);
                  padding:8px 10px; }
.chem-svg { display:block; max-width:none; flex:none; overflow:visible;
            font-family:var(--mono); fill:var(--ink); }
.chem-svg text { user-select:text; -webkit-user-select:text; }
.chem-svg.long { font-size:11.8px; }
.chem-svg.small { font-size:13.4px; }
.chem-p5{fill:var(--c-p5)} .chem-p7{fill:var(--c-p7)}
.chem-s5{fill:var(--c-s5)} .chem-s7{fill:var(--c-s7)} .chem-me{fill:var(--c-me)}
.chem-t7{fill:var(--c-t7)} .chem-cbc{fill:var(--c-cbc)} .chem-umi{fill:var(--c-umi)}
.chem-r1{fill:var(--c-r1)} .chem-r2{fill:var(--c-r2)} .chem-r3{fill:var(--c-r3)}
.chem-tso{fill:var(--c-tso)} .chem-w1{fill:var(--c-w1)}
.chem-inferred { text-decoration-line:underline; text-decoration-style:dotted;
                 text-decoration-color:var(--inf); text-underline-offset:2px; }

@media print {
  .diagram-scroll { overflow:visible; width:100%; max-width:100%; left:0; transform:none; }
  .chem-svg { max-width:100%; height:auto; }
}

.legend, .caveat { padding:12px 16px; margin:1em 0; max-width:88ch; font-size:.95rem; }
.legend { background:var(--key-bg); border-left:3px solid var(--key-rule); }
.caveat { background:var(--note-bg); border-left:3px solid var(--note-rule); }

.tw { overflow-x:auto; }
table { border-collapse:collapse; margin:.8em 0; font-size:.92rem;
        font-variant-numeric:tabular-nums; }
th, td { border:1px solid var(--rule); padding:5px 12px; text-align:left; }
th { background:var(--surface); font-weight:600; }

/* region colours: what a stretch of sequence IS */
p5{color:var(--c-p5)} p7{color:var(--c-p7)} s5{color:var(--c-s5)} s7{color:var(--c-s7)}
me{color:var(--c-me)} t7{color:var(--c-t7)} cbc{color:var(--c-cbc)} umi{color:var(--c-umi)}
r1{color:var(--c-r1)} r2{color:var(--c-r2)} r3{color:var(--c-r3)} tso{color:var(--c-tso)}
w1{color:var(--c-w1)}
/* and how well it is KNOWN -- orthogonal, nests freely inside the above */
inf { border-bottom:1px dotted var(--inf); }
</style>"""


MD_STYLE = """<style>
/* ---- extras for rendered Markdown notes --------------------------------- */
/* page.py's STYLE covers the hand-built protocol pages; these are the block
   elements only the Markdown notes use. Same tokens, so both look like one
   publication. */
.md { max-width:90ch; }
.md h4 { font-size:.95rem; font-weight:600; margin:1.6em 0 .2em; color:var(--ink); }
.md h5, .md h6 { font-size:.9rem; font-weight:600; margin:1.4em 0 .2em;
                 color:var(--ink-muted); }
.md p { max-width:82ch; }
.md ul, .md ol { max-width:82ch; padding-left:1.4em; margin:.5em 0; }
.md li { margin:.25em 0; }
.md li > ul, .md li > ol { margin:.25em 0; }
.md blockquote { margin:1em 0; padding:10px 16px; max-width:82ch;
                 background:var(--note-bg); border-left:3px solid var(--note-rule); }
.md blockquote > :first-child { margin-top:0; }
.md blockquote > :last-child { margin-bottom:0; }
.md hr { border:0; border-top:1px solid var(--rule); margin:2.4em 0; }
.md pre { max-width:100%; }
.md pre code { background:none; padding:0; font-size:.82rem; line-height:1.45; }
.md table { max-width:100%; }
.md th, .md td { vertical-align:top; }
.md td code { white-space:nowrap; }

/* breadcrumb, table of contents, footer */
.crumb { font-size:.86rem; color:var(--ink-muted); margin-bottom:1.4em; }
.crumb a { color:var(--ink-muted); }
.toc { background:var(--key-bg); border-left:3px solid var(--key-rule);
       padding:12px 18px; margin:1.4em 0 2em; max-width:82ch; font-size:.9rem; }
.toc b { display:block; font-size:.76rem; letter-spacing:.08em; text-transform:uppercase;
         color:var(--ink-muted); margin-bottom:.5em; font-weight:600; }
.toc ul { list-style:none; padding-left:0; margin:0; }
.toc li { margin:.18em 0; }
.toc li.l3 { padding-left:1.3em; font-size:.95em; }
.toc li.l3::before { content:'\\2014\\2002'; color:var(--ink-muted); }
.docfoot { margin-top:4em; padding-top:1em; border-top:1px solid var(--rule);
           font-size:.82rem; color:var(--ink-muted); }
.docfoot code { font-size:.95em; }

/* the notes' own index */
.doclist { display:grid; gap:2px; max-width:82ch; margin:1em 0; }
.doclist a { display:flex; justify-content:space-between; gap:1.5em;
             padding:7px 12px; background:var(--surface); text-decoration:none; }
.doclist a:hover { background:var(--surface-2); }
.doclist .n { color:var(--link); }
.doclist .m { color:var(--ink-muted); font-size:.85rem;
              font-variant-numeric:tabular-nums; white-space:nowrap; }
</style>"""


def head(title: str) -> str:
    """<title> plus fonts and the shared stylesheet."""
    return f"<title>{html.escape(title)}</title>\n{STYLE}"


def legend(body_html: str) -> str:
    return f'<div class="legend">\n{body_html}\n</div>'


def caveat(body_html: str) -> str:
    return f'<div class="caveat">\n{body_html}\n</div>'


def info(body_html: str) -> str:
    return f"<p><info>{body_html}</info></p>"


def table(headers, rows, scroll: bool = True) -> str:
    """A simple table; wrapped so it can scroll rather than widen the page."""
    head_html = "".join(f"<th>{h}</th>" for h in headers)
    body_html = "".join(
        "<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    t = f"<table>\n<tr>{head_html}</tr>\n{body_html}\n</table>"
    return f'<div class="tw">\n{t}\n</div>' if scroll else t
