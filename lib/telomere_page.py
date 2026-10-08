"""Small public-page renderer shared by telomere sequencing protocols."""
from __future__ import annotations

from chemdraw import Row, panel
from page import head, info


def render_page(*, title: str, notes: str, citation_html: str, panels: list[tuple[str, object, str]],
                entry: str, caveat: str | None = None) -> str:
    body = [head(title), '<div class="wrap">', f'<h1>{title}</h1>',
            f'<p class="research-notes"><a href="{notes}">Research notes</a></p>',
            info(citation_html)]
    if caveat:
        body.append(f'<div class="caveat">{caveat}</div>')
    body.append('<h2>Library construction</h2>')
    for i, (heading, drawing, caption) in enumerate(panels, 1):
        rows = drawing.rows() if hasattr(drawing, "rows") else drawing
        body.extend([f'<h3>({i}) {heading}</h3>', panel(rows, cls="long", caption=caption)])
    body.extend(['<h2>Sequencing entry</h2>',
                 panel([Row(chunks=[(entry, "r2", False)])], cls="long",
                       caption="Primer or motor entry is shown explicitly for the completed library."),
                 '</div>'])
    return "\n".join(body)
