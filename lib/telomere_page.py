"""Small public-page renderer shared by telomere sequencing protocols."""
from __future__ import annotations

from chemdraw import (Row, panel, protocol_graph_panels, workflow_from_sections,
                      workflow_panel)
from page import head, info
import seqprimers as sp


def render_page(*, title: str, notes: str, citation_html: str, panels: list[tuple[str, object, str]],
                entry: str, caveat: str | None = None, graph=None) -> str:
    if graph is not None:
        graph.require_publishable()
    wrap = ('<div class="wrap" data-protocol-complete="true">'
            if graph is not None else '<div class="wrap">')
    body = [head(title), wrap, f'<h1>{title}</h1>',
            f'<p class="research-notes"><a href="{notes}">Research notes</a></p>',
            info(citation_html)]
    if caveat:
        body.append(f'<div class="caveat">{caveat}</div>')
    body.append('<h2>Library construction</h2>')
    if graph is not None:
        body.append(protocol_graph_panels(graph, cls="long"))
    sections = []
    for heading, drawing, caption in panels:
        rows = drawing.rows() if hasattr(drawing, "rows") else drawing
        sections.append((heading, rows, caption))
    if sections and graph is None:
        body.append(workflow_panel(workflow_from_sections(sections), cls="long"))
    body.extend(['<h2>Sequencing entry</h2>',
                 panel([Row(chunks=[(entry, "r2", False)])], cls="long",
                       caption="Primer or motor entry is shown explicitly for the completed library."),
                 sp.read_spec_table([
                     sp.read_spec("Sequencing output", None, entry, "include")
                 ], "This protocol's source describes a platform-specific sequencing entry "
                    "rather than an Illumina R1/I1/I2/R2 cycle layout."),
                 '</div>'])
    return "\n".join(body)
