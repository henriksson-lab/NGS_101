"""Diagram-first page renderer for protocols producing several sequencing libraries."""
from __future__ import annotations

import seqprimers as sp
from chemdraw import Scene, annotation_rows, panel
from page import head, info


def render(m) -> str:
    parts = [head(m.TITLE), '<div class="wrap">', f'<h1>{m.TITLE}</h1>',
             f'<p class="research-notes"><a href="{m.NOTES}">Research notes</a></p>',
             info(m.SOURCE), info(m.SUMMARY)]
    if getattr(m, "CAVEAT", ""):
        parts.append(f'<div class="caveat">{m.CAVEAT}</div>')
    for heading, rows, caption in m.sections():
        parts.extend([f'<h2>{heading}</h2>', panel(rows, cls="long", caption=caption)])
    finals = tuple(getattr(m, "FINAL_LIBRARIES", ()))
    if not finals:
        raise ValueError(f"{m.TITLE}: no final sequencing libraries declared")
    for heading, library, primers, caption, intro in finals:
        if not primers:
            raise ValueError(f"{m.TITLE}: {heading} has no declared sequencing primers")
        parts.extend([f'<h2>{heading}</h2>',
                      panel([*Scene.duplex(list(library), label="library").rows(),
                             *annotation_rows(library)], cls="long", caption=caption),
                      sp.section(library, primers, intro=intro,
                                 required_roles=tuple(p.role for p in primers))])
    parts.append('</div>')
    return "\n".join(parts)
