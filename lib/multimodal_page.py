"""Diagram-first page renderer for protocols producing several sequencing libraries."""
from __future__ import annotations

import seqprimers as sp
from batch_page import oligo_rows
from chemdraw import (duplex_rows, panel, protocol_graph_panels,
                      workflow_from_sections, workflow_panel)
from page import head, info


def render(m) -> str:
    graph_factory = getattr(m, "protocol_graph", None)
    graph = graph_factory() if callable(graph_factory) else None
    if graph is not None:
        graph.require_publishable()
    graph_attr = ' data-protocol-complete="true"' if graph is not None else ""
    parts = [head(m.TITLE), f'<div class="wrap"{graph_attr}>', f'<h1>{m.TITLE}</h1>',
             f'<p class="research-notes"><a href="{m.NOTES}">Research notes</a></p>',
             info(m.SOURCE), info(m.SUMMARY)]
    if getattr(m, "CAVEAT", ""):
        parts.append(f'<div class="caveat">{m.CAVEAT}</div>')
    if callable(getattr(m, "oligos", None)):
        rows = oligo_rows(m.oligos())
        if rows:
            parts.extend(['<h2>Key oligos</h2>', '<seq>', *rows, '</seq>'])
    if graph is not None:
        workflow = None
    elif callable(getattr(m, "workflow", None)):
        workflow = m.workflow()
    else:
        sections = list(m.sections())
        workflow = workflow_from_sections(sections) if sections else None
    if graph is not None:
        parts.extend(['<h2>Reaction workflow</h2>', protocol_graph_panels(graph, cls="long")])
    elif workflow is not None:
        parts.extend(['<h2>Reaction workflow</h2>',
                      workflow_panel(workflow, cls="long")])
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
                                 read_lengths=getattr(m, "READ_LENGTHS", {}).get(heading),
                                 read_uses=getattr(m, "READ_USES", {}).get(heading),
                                 show_read_content=getattr(m, "SHOW_READ_CONTENT", True))])
    parts.append('</div>')
    return "\n".join(parts)
