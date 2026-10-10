"""Diagram-first page renderer for protocols producing several sequencing libraries."""
from __future__ import annotations

import seqprimers as sp
from chemdraw import duplex_rows, panel, workflow_from_sections, workflow_panel
from page import head, info


def render(m) -> str:
    parts = [head(m.TITLE), '<div class="wrap">', f'<h1>{m.TITLE}</h1>',
             f'<p class="research-notes"><a href="{m.NOTES}">Research notes</a></p>',
             info(m.SOURCE), info(m.SUMMARY)]
    if getattr(m, "CAVEAT", ""):
        parts.append(f'<div class="caveat">{m.CAVEAT}</div>')
    if callable(getattr(m, "oligos", None)):
        rows = list(m.oligos())
        if rows:
            parts.extend(['<h2>Key oligos</h2>', '<seq>', *rows, '</seq>'])
    sections = list(m.sections())
    if sections:
        parts.extend(['<h2>Reaction workflow</h2>',
                      workflow_panel(workflow_from_sections(sections), cls="long")])
    finals = tuple(getattr(m, "FINAL_LIBRARIES", ()))
    if not finals:
        raise ValueError(f"{m.TITLE}: no final sequencing libraries declared")
    for heading, library, primers, caption, intro in finals:
        if not primers:
            raise ValueError(f"{m.TITLE}: {heading} has no declared sequencing primers")
        parts.extend([f'<h2>{heading}</h2>',
                      panel(duplex_rows(library, label="library"),
                            cls="long", caption=caption),
                      sp.section(library, primers, intro=intro,
                                 required_roles=tuple(p.role for p in primers),
                                 read_lengths=getattr(m, "READ_LENGTHS", {}).get(heading))])
    parts.append('</div>')
    return "\n".join(parts)
