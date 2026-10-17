"""Small public-page renderer shared by telomere sequencing protocols."""
from __future__ import annotations

from chemdraw import (Row, panel, protocol_graph_panels, workflow_from_sections,
                      workflow_panel)
from page import head, info
import seqprimers as sp
from batch_page import oligo_rows
from protocol_graph import AnnotationSlot
from protocol_views import render_annotations


def render_page(*, title: str, notes: str, citation_html: str, panels: list[tuple[str, object, str]],
                entry: str, caveat: str | None = None, graph=None, oligos=()) -> str:
    if graph is not None:
        graph.require_publishable()
    wrap = ('<div class="wrap" data-protocol-complete="true">'
            if graph is not None else '<div class="wrap">')
    body = [head(title), wrap, f'<h1>{title}</h1>',
            f'<p class="research-notes"><a href="{notes}">Research notes</a></p>',
            info(citation_html)]
    def annotations(slot):
        markup = render_annotations(graph, slot) if graph is not None else ""
        if markup:
            body.append(markup)
    if caveat:
        body.append(f'<div class="caveat">{caveat}</div>')
    annotations(AnnotationSlot.INTRO)
    rows = oligo_rows(oligos)
    if rows:
        body.extend(['<h2>Key oligos</h2>', '<seq>', *rows, '</seq>'])
    body.append('<h2>Library construction</h2>')
    if graph is not None:
        annotations(AnnotationSlot.BEFORE_WORKFLOW)
        body.append(protocol_graph_panels(graph, cls="long"))
        annotations(AnnotationSlot.AFTER_WORKFLOW)
    sections = []
    for heading, drawing, caption in panels:
        rows = drawing.rows() if hasattr(drawing, "rows") else drawing
        sections.append((heading, rows, caption))
    if sections and graph is None:
        body.append(workflow_panel(workflow_from_sections(sections), cls="long"))
    body.append('<h2>Sequencing entry</h2>')
    annotations(AnnotationSlot.SEQUENCING)
    body.extend([panel([Row(chunks=[(entry, "r2", False)])], cls="long",
                       caption="Primer or motor entry is shown explicitly for the completed library."),
                 sp.read_spec_table([
                     sp.read_spec("Sequencing output", None, entry, "include")
                 ], "This protocol's source describes a platform-specific sequencing entry "
                    "rather than an Illumina R1/I1/I2/R2 cycle layout."),
                 '</div>'])
    return "\n".join(body)
