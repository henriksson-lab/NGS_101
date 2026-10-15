"""Diagram-first renderer shared by compact, source-backed protocol modules."""
from __future__ import annotations

import html
import re

import seqprimers as sp
from chemdraw import (duplex_rows, panel, protocol_graph_panels,
                      workflow_from_sections, workflow_panel)
from page import head, info
from protocol_graph import AnnotationSlot
from protocol_views import render_annotations


def oligo_rows(rows):
    """Normalize ordering rows and reject a strand written backwards.

    Source diagrams sometimes draw a paired strand left-to-right as 3′→5′.  An ordering
    list must instead preserve the oligo's written 5′→3′ sequence; accepting the former
    here would make an attractive but dangerous copy target.
    """
    out = []
    reverse_order = re.compile(
        r"3[′']\s*[-–]\s*(?:/[A-Za-z0-9]+/|\[[^]]+\]|[ACGTUNRYSWKMBDHVX*+\s]){8,}"
        r"\s*[-–]\s*5[′']", re.I)
    for row in rows:
        plain = html.unescape(re.sub(r"<[^>]+>", "", row))
        if reverse_order.search(plain):
            raise ValueError("ordering oligo is written 3-prime to 5-prime: " + plain)
        out.append(row if row.lstrip().startswith("<p") else f"<p>{row}</p>")
    return out


def render(m) -> str:
    graph_factory = getattr(m, "protocol_graph", None)
    graph = graph_factory() if callable(graph_factory) else None
    if graph is not None:
        graph.require_publishable()
    parts = [head(m.TITLE), '<div class="wrap">', f'<h1>{m.TITLE}</h1>',
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
        workflow = (workflow_from_sections(
            sections, initial_rows=getattr(m, "INITIAL_ROWS", None),
            initial_name=getattr(m, "INITIAL_NAME", "Starting material"))
            if sections else None)
    if graph is not None:
        annotations(AnnotationSlot.BEFORE_WORKFLOW)
        parts[1] = '<div class="wrap" data-protocol-complete="true">'
        parts.extend(['<h2>Reaction workflow</h2>', protocol_graph_panels(graph, cls="long")])
        annotations(AnnotationSlot.AFTER_WORKFLOW)
    elif workflow is not None:
        parts.extend(['<h2>Reaction workflow</h2>', workflow_panel(workflow, cls="long")])
    for title, markup in getattr(m, "TOPOLOGY_PANELS", ()):
        parts.extend([f'<h2>{title}</h2>', markup])
    plan = (next(iter(graph.sequencing.values()))
            if graph is not None and len(graph.sequencing) == 1 else None)
    # A declared None means this page draws a topology (for example a circular
    # SMRTbell) instead of a linear final-library duplex.  A graph sequencing
    # circuit is a read path, not necessarily that physical final molecule.
    topology_endpoint = (hasattr(m, "FINAL_LIBRARY") and m.FINAL_LIBRARY is None
                         and not hasattr(m, "FINAL_CAPTION"))
    lib = (None if topology_endpoint else
           ((plan.library if plan is not None else None)
            or getattr(m, "FINAL_LIBRARY", None)))
    primers = ((plan.primers if plan is not None else ())
               or getattr(m, "SEQ_PRIMERS", ()))
    if lib is not None:
        # Bracketed role tokens deliberately stand for unavailable molecular sequence.
        # Scene still owns placement; excluding those prose columns from base-pair
        # validation avoids pretending that a literal "[adapter]" has a complement.
        unknown = tuple(s.name + "'" for s in lib if s.is_role_token())
        final_rows = (m.final_rows() if callable(getattr(m, "final_rows", None))
                      else duplex_rows(lib, label="library", unpaired=unknown))
        parts.extend(['<h2>Final sequencing library</h2>',
                      panel(final_rows, cls="long",
                            caption=m.FINAL_CAPTION)])
    for title, state, caption in getattr(m, "EXTRA_MOLECULAR_STATES", ()):
        unknown = tuple(s.name + "'" for s in state if s.is_role_token())
        parts.extend([f'<h2>{title}</h2>',
                      panel(duplex_rows(state, label="library", unpaired=unknown),
                            cls="long", caption=caption)])
    if primers:
        annotations(AnnotationSlot.SEQUENCING)
        acquired_roles = ((plan.acquired_roles if plan is not None else ())
                          or getattr(m, "ACQUIRED_ROLES", tuple(p.role for p in primers)))
        show_read_content = getattr(m, "SHOW_READ_CONTENT", True)
        parts.append(sp.section(lib, primers, intro=m.SEQUENCING_INTRO,
                                required_roles=acquired_roles,
                                read_lengths=(plan.read_lengths if plan is not None and
                                              plan.read_lengths else getattr(m, "READ_LENGTHS", None)),
                                read_uses=(plan.read_uses if plan is not None and
                                           plan.read_uses else getattr(m, "READ_USES", None)),
                                show_read_content=show_read_content))
        if not show_read_content:
            specs = getattr(m, "READ_SPECS", ())
            if not specs:
                raise ValueError(
                    f"{m.TITLE}: suppressed geometry-derived reads need READ_SPECS")
            parts.append(sp.read_spec_table(specs))
        unavailable = getattr(m, "SEQUENCING_UNAVAILABLE", "")
        if unavailable:
            parts.append(sp.unavailable_diagram(lib, unavailable))
    else:
        annotations(AnnotationSlot.SEQUENCING)
        ending = getattr(m, "SEQUENCING_ENDING", "")
        if not ending:
            raise ValueError(f"{m.TITLE}: page has neither sequencing primers nor a declared endpoint")
        parts.extend(['<h2>Sequencing entry</h2>', info(ending)])
        unavailable = getattr(m, "SEQUENCING_UNAVAILABLE", "")
        if unavailable:
            parts.append(sp.unavailable_diagram(lib, unavailable))
        specs = getattr(m, "READ_SPECS", ())
        if specs:
            parts.append(sp.read_spec_table(specs))
        elif getattr(m, "SEQUENCING_MECHANISM", "") == "nanopore_motor":
            parts.append(sp.nanopore_read_content_table(
                lib, getattr(m, "MOTOR_ADAPTER_NAME", "motor adapter")))
        else:
            parts.append(sp.unresolved_read_content_table(
                lib, unavailable or "the exact sequencing-primer site is not represented",
                read_lengths=getattr(m, "READ_LENGTHS", None),
                acquired_roles=getattr(m, "ACQUIRED_ROLES", None)))
    parts.extend(['</div>'])
    return "\n".join(parts)
