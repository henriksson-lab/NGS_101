"""Diagram-first renderer shared by compact, source-backed protocol modules."""
from __future__ import annotations

import html
import re

import seqprimers as sp
from chemdraw import (Construct, duplex_rows, panel, protocol_graph_panels,
                      workflow_from_sections, workflow_panel)
from page import head, info
from protocol_graph import AnnotationSlot
from protocol_views import render_annotations, sequencing_limitation


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
    plans = _distinct_plans(graph)
    if len(plans) > 1 and any(_exact(plan) for _, plan in plans):
        _render_terminal_sequencing(m, graph, plans, parts, annotations)
    else:
        _render_page_sequencing(m, graph, plans, parts, annotations)
    parts.extend(['</div>'])
    return "\n".join(parts)


def _distinct_plans(graph):
    """``(material ids, plan)`` per distinct graph sequencing plan, in graph order."""
    plans: list[tuple[list[str], object]] = []
    if graph is None:
        return plans
    for material_id, plan in graph.sequencing.items():
        for ids, known in plans:
            if known is plan or known == plan:
                ids.append(material_id)
                break
        else:
            plans.append(([material_id], plan))
    return plans


def _plain(text: str) -> str:
    return " ".join(html.unescape(text).split())


def _exact(plan) -> bool:
    return plan.library is not None and bool(plan.primers)


def _unknown_columns(state):
    # Bracketed role tokens deliberately stand for unavailable molecular sequence.
    # Scene still owns placement; excluding those prose columns from base-pair
    # validation avoids pretending that a literal "[adapter]" has a complement.
    return tuple(s.name + "'" for s in state if s.is_role_token())


def _read_specs_note(m, has_geometry: bool) -> str:
    note = getattr(m, "READ_SPECS_NOTE", "")
    if note:
        return note
    return sp.PARTIAL_GEOMETRY_NOTE if has_geometry else ""


def _extra_states(m, parts, drawn=()):
    for title, state, caption in getattr(m, "EXTRA_MOLECULAR_STATES", ()):
        if any(state is library for library in drawn):
            continue  # already drawn from the graph's own sequencing plan
        parts.extend([f'<h2>{title}</h2>',
                      panel(duplex_rows(state, label="library",
                                        unpaired=_unknown_columns(state)),
                            cls="long", caption=caption)])


def _render_terminal_sequencing(m, graph, plans, parts, annotations):
    """One final-library and sequencing section per distinct graph plan.

    Incomplete branches retain their own limitation; a neighbouring exact branch must
    never supply their sequencing geometry or make them disappear into a page fallback.
    """
    if not getattr(m, "SHOW_READ_CONTENT", True):
        raise ValueError(f"{m.TITLE}: suppressed read content needs one sequencing plan")
    captions = getattr(m, "FINAL_CAPTIONS", {})
    final = getattr(m, "FINAL_LIBRARY", None)
    annotations(AnnotationSlot.SEQUENCING)
    if getattr(m, "SEQUENCING_INTRO", ""):
        parts.append(info(m.SEQUENCING_INTRO))
    for ids, plan in plans:
        label = "; ".join(graph.materials[i].label for i in ids)
        if not _exact(plan):
            parts.append(f'<h2>Final sequencing library: {html.escape(label)}</h2>')
            state = plan.library if plan.library is not None else graph.materials[ids[0]].state
            if isinstance(state, Construct):
                caption = ("INFERRED — " if any(s.inferred for s in state) else "") + label
                parts.append(panel(duplex_rows(state, label="library",
                                               unpaired=_unknown_columns(state)),
                                   cls="long", caption=caption))
            elif callable(getattr(state, "rows", None)):
                parts.append(panel(state.rows(), cls="long", caption=label))
            parts.append(sequencing_limitation(plan))
            parts.append(sp.unresolved_read_content_table(
                state if isinstance(state, Construct) else None,
                plan.limitation, read_lengths=plan.read_lengths or None,
                acquired_roles=plan.acquired_roles or None))
            continue
        caption = next((captions[i] for i in ids if i in captions), None)
        if caption is None:
            caption = next((extra_caption for _, state, extra_caption
                            in getattr(m, "EXTRA_MOLECULAR_STATES", ())
                            if state is plan.library), None)
        if caption is None:
            caption = (getattr(m, "FINAL_CAPTION", "") if plan.library is final
                       else "") or label
        parts.extend([f'<h2>Final sequencing library: {html.escape(label)}</h2>',
                      panel(duplex_rows(plan.library, label="library",
                                        unpaired=_unknown_columns(plan.library)),
                            cls="long", caption=caption)])
        parts.append(sp.section(plan.library, plan.primers,
                                heading=f"Sequencing primers: {label}",
                                required_roles=plan.required_roles,
                                read_lengths=plan.read_lengths or None,
                                read_uses=plan.read_uses or None,
                                unplaced_reads=plan.unplaced_reads,
                                acquisition_unspecified=plan.acquisition_unspecified))
        limitation = sequencing_limitation(plan)
        if limitation:
            parts.append(limitation)
    _extra_states(m, parts, drawn=tuple(plan.library for _, plan in plans))


def _render_page_sequencing(m, graph, plans, parts, annotations):
    """The page-level sequencing section; a single distinct graph plan overrides it."""
    plan = plans[0][1] if len(plans) == 1 else None
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
        final_rows = (m.final_rows() if callable(getattr(m, "final_rows", None))
                      else duplex_rows(lib, label="library", unpaired=_unknown_columns(lib)))
        parts.extend(['<h2>Final sequencing library</h2>',
                      panel(final_rows, cls="long",
                            caption=m.FINAL_CAPTION)])
    _extra_states(m, parts)
    unavailable = getattr(m, "SEQUENCING_UNAVAILABLE", "")
    # Limitations recorded on graph plans are shown unless the page already states them.
    stated = _plain("\n".join(parts) + "\n" + "\n".join(
        str(getattr(m, name, "") or "") for name in (
            "SEQUENCING_UNAVAILABLE", "SEQUENCING_ENDING", "SEQUENCING_INTRO",
            "FINAL_CAPTION", "READ_SPECS_NOTE"))
        + "\n" + "\n".join(getattr(row, "content", "")
                             for row in getattr(m, "READ_SPECS", ())))
    limitations = "".join(
        sequencing_limitation(p, "" if len(plans) == 1 else
                              "; ".join(graph.materials[i].label for i in ids))
        for ids, p in plans
        if p.limitation.strip() and _plain(p.limitation) not in stated)
    if primers:
        annotations(AnnotationSlot.SEQUENCING)
        exact_plan = plan is not None and plan.primers
        acquisition_unspecified = bool(exact_plan and plan.acquisition_unspecified)
        unplaced = plan.unplaced_reads if exact_plan else ()
        acquired_roles = (() if acquisition_unspecified else
                          (plan.acquired_roles if plan is not None else ())
                          or getattr(m, "ACQUIRED_ROLES", tuple(p.role for p in primers)))
        show_read_content = getattr(m, "SHOW_READ_CONTENT", True)
        parts.append(sp.section(lib, primers, intro=m.SEQUENCING_INTRO,
                                required_roles=acquired_roles,
                                read_lengths=(plan.read_lengths if plan is not None and
                                              plan.read_lengths else getattr(m, "READ_LENGTHS", None)),
                                read_uses=(plan.read_uses if plan is not None and
                                           plan.read_uses else getattr(m, "READ_USES", None)),
                                show_read_content=show_read_content,
                                unplaced_reads=unplaced,
                                acquisition_unspecified=acquisition_unspecified))
        if not show_read_content:
            specs = getattr(m, "READ_SPECS", ())
            if not specs:
                raise ValueError(
                    f"{m.TITLE}: suppressed geometry-derived reads need READ_SPECS")
            parts.append(sp.read_spec_table(specs, _read_specs_note(m, True)))
        if limitations:
            parts.append(limitations)
        if unavailable:
            parts.append(sp.unavailable_diagram(lib, unavailable))
    else:
        annotations(AnnotationSlot.SEQUENCING)
        ending = getattr(m, "SEQUENCING_ENDING", "")
        if not ending:
            raise ValueError(f"{m.TITLE}: page has neither sequencing primers nor a declared endpoint")
        parts.extend(['<h2>Sequencing entry</h2>', info(ending)])
        if limitations:
            parts.append(limitations)
        if unavailable:
            parts.append(sp.unavailable_diagram(lib, unavailable))
        specs = getattr(m, "READ_SPECS", ())
        if specs:
            parts.append(sp.read_spec_table(specs, _read_specs_note(m, False)))
        elif getattr(m, "SEQUENCING_MECHANISM", "") == "nanopore_motor":
            parts.append(sp.nanopore_read_content_table(
                lib, getattr(m, "MOTOR_ADAPTER_NAME", "motor adapter")))
        else:
            parts.append(sp.unresolved_read_content_table(
                lib, unavailable or "the exact sequencing-primer site is not represented",
                read_lengths=getattr(m, "READ_LENGTHS", None),
                acquired_roles=getattr(m, "ACQUIRED_ROLES", None)))
