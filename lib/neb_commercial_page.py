"""Diagram-first renderer for NEB commercial protocol modules."""
from __future__ import annotations

import seqprimers as sp
from chemdraw import (Construct, Scene, annotation_rows, duplex_rows, panel,
                      protocol_graph_panels, strand_row, MolecularState, Row,
                      Workflow, workflow_panel)
from page import head, info
from protocol_graph import AnnotationSlot
from protocol_views import render_annotations


def _rows(obj):
    if isinstance(obj, Scene):
        return obj.rows()
    if isinstance(obj, Construct):
        return [strand_row(obj), *annotation_rows(obj)]
    return obj


def render(module) -> str:
    graph_factory = getattr(module, "protocol_graph", None)
    graph = graph_factory() if callable(graph_factory) else None
    if graph is not None:
        graph.require_publishable()
    chunks = [head(module.TITLE), '<div class="wrap">', f'<h1>{module.TITLE}</h1>',
              f'<p class="research-notes"><a href="{module.NOTES}">Research notes</a></p>',
              info(module.SOURCE), info(module.SUMMARY)]
    caveat = getattr(module, "CAVEAT", "")
    if caveat:
        chunks.append(f'<div class="caveat"><b>INFERRED boundary.</b> {caveat}</div>')
    if graph is not None:
        chunks.append(render_annotations(graph, AnnotationSlot.INTRO))
    if graph is not None:
        chunks[1] = '<div class="wrap" data-protocol-complete="true">'
        chunks.append(render_annotations(graph, AnnotationSlot.BEFORE_WORKFLOW))
        chunks.extend(['<h2>Reaction workflow</h2>',
                       protocol_graph_panels(graph, cls="long")])
        chunks.append(render_annotations(graph, AnnotationSlot.AFTER_WORKFLOW))
    else:
        workflow = Workflow(MolecularState(
            getattr(module, "INITIAL_NAME", "Starting material"),
            tuple(getattr(module, "INITIAL_ROWS", (Row(chunks=[("[starting material]", None, False)]),)))))
        for step in module.STEPS:
            heading, obj, caption, *geometry = step
            rows = _rows(obj)
            if geometry:
                workflow.ligate_adapter(heading, geometry[0], rows, name=heading, note=caption)
            else:
                workflow.react(heading, rows, name=heading, note=caption)
        if workflow.reactions:
            chunks.extend(['<h2>Reaction workflow</h2>',
                           workflow_panel(workflow, cls="small")])
    lib = module.final_library()
    chunks.extend(['<h2>Final sequencing molecule</h2>',
                   panel(duplex_rows(lib, label="library"),
                         cls="long", caption=module.FINAL_CAPTION),
                   sp.section(lib, module.SEQ_PRIMERS, intro=module.READOUT,
                              required_roles=module.REQUIRED_ROLES,
                              read_lengths=getattr(module, "READ_LENGTHS", None))])
    if graph is not None:
        chunks.append(render_annotations(graph, AnnotationSlot.SEQUENCING))
    chunks.append('</div>')
    return "\n".join(chunks)
