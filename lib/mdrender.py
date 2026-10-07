"""
Dependency-free Markdown -> HTML for this repo's .md notes.

The repo has no third-party dependencies and none of python-markdown, mistune,
commonmark or markdown-it is installed, so this is a small renderer covering
exactly the subset the notes actually use. Surveyed across 34 files / 7,288
lines: ATX headings, GFM tables with alignment rows, fenced code, blockquotes,
nested ordered/unordered lists, thematic breaks, and inline
code/bold/italic/links/autolinks/bare URLs.

Two deliberate departures from CommonMark, each forced by what is in the notes:

1. **No indented code blocks.** Every code block in the notes is fenced (156
   fences). All 84 four-space indents are list continuations, so an indented
   line is treated as continuation text and never as code.

2. **Raw HTML is escaped, except a small inline allowlist.** The notes use
   angle brackets as *metavariables* far more often than as markup --
   `<cbc>`, `<seqid>`, `<align>`, `<inf>`, `<name>`, `<id>` -- and those must
   render literally rather than vanish into the DOM. Only `<br> <sub> <sup>
   <i> <b> <em> <strong>` pass through, which covers every genuine use
   (line breaks inside table cells, T<sub>m</sub>, author footnote markers).

Links to `.md` files are rewritten to `.html` so the generated site navigates
the same way the sources do.
"""

from __future__ import annotations

import html as _html
import re

ALLOWED_INLINE_HTML = ("br", "sub", "sup", "i", "b", "em", "strong")

_MARKER = re.compile(r"^(\s*)(?:([-*+])|(\d+)[.)])\s+")
_FENCE = re.compile(r"^(\s*)(`{3,}|~{3,})\s*([\w+#-]*)\s*$")
_HEADING = re.compile(r"^(#{1,6})\s+(.*?)(?:\s+#+)?\s*$")
_HR = re.compile(r"^\s{0,3}([-*_])\s*(?:\1\s*){2,}$")
_ALLOWED_RE = re.compile(r"&lt;(/?)(" + "|".join(ALLOWED_INLINE_HTML) + r")\s*/?&gt;")


def slug(text: str) -> str:
    """A heading's anchor id: markup and emoji stripped, spaces to hyphens."""
    s = re.sub(r"`|\*+|_", "", text)
    s = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", s)          # links -> their label
    s = re.sub(r"[^\w\s-]", "", s, flags=re.UNICODE)        # drops emoji and punctuation
    return re.sub(r"[\s_]+", "-", s.strip().lower()) or "section"


def rewrite_link(url: str, md2html: bool = True) -> str:
    """Point .md links at their generated .html sibling, anchors preserved."""
    if not md2html or url.startswith(("http://", "https://", "mailto:", "#")):
        return url
    return re.sub(r"\.md(?=$|#)", ".html", url)


class _Store:
    """Holds finished HTML behind placeholders so later passes cannot touch it."""

    def __init__(self) -> None:
        self.items: list[str] = []

    def put(self, html_str: str) -> str:
        self.items.append(html_str)
        return f"\x00{len(self.items) - 1}\x01"

    def restore(self, text: str) -> str:
        # Stored HTML may itself contain placeholders (a link label holding a
        # code span), so expand until none are left.
        for _ in range(12):
            new = re.sub(r"\x00(\d+)\x01", lambda m: self.items[int(m.group(1))], text)
            if new == text:
                return new
            text = new
        return text


def _bare_url(m: re.Match, store: _Store) -> str:
    url = m.group(0).rstrip(".,;:!?)")
    tail = m.group(0)[len(url):]
    return store.put(f'<a href="{url}">{url}</a>') + tail


def inline(text: str, store: _Store, md2html: bool = True) -> str:
    """Inline markup, in an order chosen so each pass cannot corrupt the next."""
    # 1. code spans first: nothing inside a backtick run is markup
    text = re.sub(r"(`+)(.+?)\1",
                  lambda m: store.put(f"<code>{_html.escape(m.group(2), quote=False)}</code>"),
                  text, flags=re.S)
    # 2. everything that is left is literal text
    text = _html.escape(text, quote=False)
    # 3. let the allowlisted inline tags back through
    text = _ALLOWED_RE.sub(lambda m: store.put(f"<{m.group(1)}{m.group(2)}>"), text)
    # 4. <https://...> autolinks (now escaped, so match the escaped form)
    text = re.sub(r"&lt;(https?://[^\s&]+)&gt;",
                  lambda m: store.put(f'<a href="{m.group(1)}">{m.group(1)}</a>'), text)
    # 5. emphasis, strongest marker first
    text = re.sub(r"\*\*(?=\S)(.+?)(?<=\S)\*\*", r"<strong>\1</strong>", text, flags=re.S)
    text = re.sub(r"(?<![\w*])\*(?=[^\s*])([^*]+?)(?<=\S)\*(?![\w*])", r"<em>\1</em>", text)
    text = re.sub(r"(?<![\w_])_(?=\S)([^_]+?)(?<=\S)_(?![\w_])", r"<em>\1</em>", text)
    # 6. [label](url) -- the label has already been through 1-5
    text = re.sub(
        r"\[([^\]]*)\]\(([^)\s]+)\)",
        lambda m: store.put(
            f'<a href="{_html.escape(rewrite_link(m.group(2), md2html), quote=True)}">'
            f"{m.group(1)}</a>"),
        text)
    # 7. remaining bare URLs
    text = re.sub(r"(?<![\w@.\"'/=])https?://[^\s<>()\[\]\"']+",
                  lambda m: _bare_url(m, store), text)
    return text


def _is_align_row(line: str, ncols: int) -> bool:
    cells = _split_row(line)
    return (len(cells) == ncols and bool(cells)
            and all(re.fullmatch(r":?-+:?", c) for c in cells))


def _split_row(line: str) -> list[str]:
    s = line.strip()
    if s.startswith("|"):
        s = s[1:]
    if s.endswith("|") and not s.endswith("\\|"):
        s = s[:-1]
    return [c.strip() for c in s.split("|")]


def _align_of(cell: str) -> str:
    if cell.startswith(":") and cell.endswith(":"):
        return ' style="text-align:center"'
    if cell.endswith(":"):
        return ' style="text-align:right"'
    return ""


def render(md: str, *, md2html: bool = True, heading_ids: bool = True,
           headings: list[tuple[int, str, str]] | None = None) -> str:
    """Render markdown to body HTML. `headings` collects (level, text, id)."""
    store = _Store()
    out: list[str] = []
    lines = md.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    i, n = 0, len(lines)

    while i < n:
        line = lines[i]

        if not line.strip():
            i += 1
            continue

        m = _FENCE.match(line)
        if m:
            close, lang, buf = m.group(2)[0], m.group(3), []
            i += 1
            while i < n and not re.match(rf"^\s*{re.escape(close)}{{3,}}\s*$", lines[i]):
                buf.append(lines[i])
                i += 1
            i += 1                                        # consume the closing fence
            cls = f' class="language-{lang}"' if lang else ""
            body = _html.escape("\n".join(buf), quote=False)
            out.append(f"<pre><code{cls}>{body}</code></pre>")
            continue

        if _HR.match(line):
            out.append("<hr>")
            i += 1
            continue

        m = _HEADING.match(line)
        if m:
            lvl, raw = len(m.group(1)), m.group(2)
            sid = slug(raw)
            if headings is not None:
                headings.append((lvl, raw, sid))
            attr = f' id="{sid}"' if heading_ids else ""
            out.append(f"<h{lvl}{attr}>{store.restore(inline(raw, store, md2html))}</h{lvl}>")
            i += 1
            continue

        if "|" in line and i + 1 < n and _is_align_row(lines[i + 1], len(_split_row(line))):
            hdr = _split_row(line)
            aligns = [_align_of(c) for c in _split_row(lines[i + 1])]
            i += 2
            rows = []
            while i < n and "|" in lines[i] and lines[i].strip():
                rows.append(_split_row(lines[i]))
                i += 1
            th = "".join(f"<th{a}>{store.restore(inline(c, store, md2html))}</th>"
                         for c, a in zip(hdr, aligns))
            body = []
            for r in rows:
                tds = "".join(
                    f"<td{aligns[k] if k < len(aligns) else ''}>"
                    f"{store.restore(inline(c, store, md2html))}</td>"
                    for k, c in enumerate(r))
                body.append(f"<tr>{tds}</tr>")
            out.append('<div class="tw">\n<table>\n<tr>' + th + "</tr>\n"
                       + "\n".join(body) + "\n</table>\n</div>")
            continue

        if line.lstrip().startswith(">"):
            buf = []
            while i < n and (lines[i].lstrip().startswith(">") or
                             (lines[i].strip() and buf)):
                buf.append(re.sub(r"^\s*>\s?", "", lines[i]))
                i += 1
            out.append("<blockquote>\n"
                       + render("\n".join(buf), md2html=md2html,
                                heading_ids=False, headings=headings)
                       + "\n</blockquote>")
            continue

        if _MARKER.match(line):
            html_str, i = _list(lines, i, md2html, headings)
            out.append(html_str)
            continue

        buf = []
        while i < n and lines[i].strip() and not _FENCE.match(lines[i]) \
                and not _HEADING.match(lines[i]) and not _HR.match(lines[i]) \
                and not lines[i].lstrip().startswith(">") and not _MARKER.match(lines[i]):
            buf.append(lines[i].strip())
            i += 1
        if buf:
            out.append(f"<p>{store.restore(inline(' '.join(buf), store, md2html))}</p>")

    return "\n".join(out)


def _list(lines: list[str], i: int, md2html: bool,
          headings: list | None) -> tuple[str, int]:
    """Parse one list, recursing for nested content. Returns (html, next_i)."""
    m0 = _MARKER.match(lines[i])
    assert m0 is not None
    base = len(m0.group(1))
    ordered = m0.group(3) is not None
    start = int(m0.group(3)) if ordered else 1
    items: list[list[str]] = []
    cinds: list[int] = []
    n = len(lines)

    while i < n:
        line = lines[i]
        if not line.strip():
            j = i
            while j < n and not lines[j].strip():
                j += 1
            if j >= n:
                break
            ind = len(lines[j]) - len(lines[j].lstrip())
            mm = _MARKER.match(lines[j])
            if ind > base or (mm and len(mm.group(1)) == base):
                if items:
                    items[-1].append("")
                i = j
                continue
            break
        ind = len(line) - len(line.lstrip())
        mm = _MARKER.match(line)
        if mm and len(mm.group(1)) == base:
            # a marker of the other kind at this level starts a separate list
            if items and ((mm.group(3) is not None) != ordered):
                break
            items.append([line[mm.end():]])
            cinds.append(mm.end())
            i += 1
        elif ind > base and items:
            cind = cinds[-1]
            items[-1].append(line[cind:] if ind >= cind else line.lstrip())
            i += 1
        elif items and ind <= base and not mm:
            items[-1].append(line.strip())           # lazy continuation
            i += 1
        else:
            break

    lis = []
    for content in items:
        inner = render("\n".join(content), md2html=md2html,
                       heading_ids=False, headings=headings).strip()
        one = re.fullmatch(r"<p>(.*)</p>", inner, flags=re.S)
        lis.append(f"<li>{one.group(1)}</li>" if one else f"<li>\n{inner}\n</li>")
    tag = "ol" if ordered else "ul"
    attr = f' start="{start}"' if ordered and start != 1 else ""
    return f"<{tag}{attr}>\n" + "\n".join(lis) + f"\n</{tag}>", i
