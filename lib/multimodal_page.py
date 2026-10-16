"""Diagram-first page renderer for protocols producing several sequencing libraries."""
from __future__ import annotations

import seqprimers as sp
from batch_page import oligo_rows
from chemdraw import (duplex_rows, panel, protocol_graph_panels,
                      workflow_from_sections, workflow_panel)
from page import head, info
from protocol_graph import AnnotationSlot
from protocol_views import render_annotations


def render(m) -> str:
    graph_factory = getattr(m, "protocol_graph", None)
    graph = graph_factory() if callable(graph_factory) else None
    if graph is not None:
        graph.require_publishable()
    graph_attr = ' data-protocol-complete="true"' if graph is not None else ""
    parts = [head(m.TITLE), f'<div class="wrap"{graph_attr}>', f'<h1>{m.TITLE}</h1>',
             f'<p class="research-notes"><a href="{m.NOTES}">Research notes</a></p>',
             info(m.SOURCE), info(m.SUMMARY)]
    def annotations(slot):
        if graph is not None:
            markup = render_annotations(graph, slot)
            if markup:
                parts.append(markup)
    if getattr(m, "CAVEAT", ""):
        parts.append(f'<div class="caveat">{m.CAVEAT}</div>')
    annotations(AnnotationSlot.INTRO)
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
        annotations(AnnotationSlot.BEFORE_WORKFLOW)
        parts.extend(['<h2>Reaction workflow</h2>', protocol_graph_panels(graph, cls="long")])
        annotations(AnnotationSlot.AFTER_WORKFLOW)
    elif workflow is not None:
        parts.extend(['<h2>Reaction workflow</h2>',
                      workflow_panel(workflow, cls="long")])
    finals = tuple(getattr(m, "FINAL_LIBRARIES", ()))
    if not finals:
        raise ValueError(f"{m.TITLE}: no final sequencing libraries declared")
    for index, (heading, library, primers, caption, intro) in enumerate(finals):
        if index == 0:
            annotations(AnnotationSlot.SEQUENCING)
        plan = None
        if graph is not None:
            matches = [candidate for candidate in graph.sequencing.values()
                       if candidate.library is library or candidate.library == library]
            if len(matches) > 1:
                raise ValueError(f"{m.TITLE}: multiple graph readouts match {heading}")
            plan = matches[0] if matches else None
        if plan is not None:
            library = plan.library
            primers = plan.primers
        specs = getattr(m, "READ_SPECS", {}).get(heading, ())
        limitations = getattr(m, "SEQUENCING_UNAVAILABLE", {})
        limitation = (plan.limitation if plan is not None else "") or (
            limitations.get(heading, "") if isinstance(limitations, dict) else limitations)
        if not primers and not (limitation and specs):
            raise ValueError(f"{m.TITLE}: {heading} needs sequencing primers or a limitation and source-declared reads")
        parts.extend([f'<h2>{heading}</h2>',
                      panel(duplex_rows(library, label="library"),
                            cls="long", caption=caption)])
        if primers:
            parts.append(sp.section(library, primers, intro=intro,
                                 required_roles=(plan.acquired_roles if plan is not None and
                                                 plan.acquired_roles else tuple(p.role for p in primers)),
                                 read_lengths=(plan.read_lengths if plan is not None and
                                               plan.read_lengths else
                                               getattr(m, "READ_LENGTHS", {}).get(heading)),
                                 read_uses=(plan.read_uses if plan is not None and
                                            plan.read_uses else
                                            getattr(m, "READ_USES", {}).get(heading)),
                                 show_read_content=not specs and getattr(m, "SHOW_READ_CONTENT", True)))
        if limitation:
            parts.append(sp.unavailable_diagram(library, limitation))
        if specs:
            parts.append(sp.read_spec_table(specs))
    parts.append('</div>')
    return "\n".join(parts)
