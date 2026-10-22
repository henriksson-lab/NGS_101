"""Render graph-targeted explanations; chemistry validation stays independent."""
from __future__ import annotations

import html
from dataclasses import dataclass
from typing import Iterable

from chemdraw import Construct, MolecularState, Row, Scene, Segment, oligo, panel, protocol_graph_panels
from page import caveat, info, table
import seqprimers as sp
from protocol_graph import (AnnotationSlot, DetailStyle, DetailView, Highlight, ProtocolGraph,
                            SourceDiscrepancy)


@dataclass(frozen=True)
class OligoEntry:
    name: str
    segments: tuple[Segment, ...]
    mods: str = ""


@dataclass(frozen=True)
class OligoSection:
    title: str
    entries: tuple[OligoEntry, ...]
    note_html: str = ""
    multiline: bool = False

    def render(self) -> str:
        if not self.entries:
            raise ValueError("oligo section needs at least one entry")
        sep = "\n" if self.multiline else ""
        opening = f"<h2>{html.escape(self.title)}</h2>" + ("\n<seq>\n" if self.multiline else "<seq>")
        closing = "\n</seq>\n" if self.multiline else "</seq>"
        body = sep.join(oligo(item.name, item.segments, mods=item.mods)
                        for item in self.entries)
        return opening + body + closing + (info(self.note_html) if self.note_html else "")


@dataclass(frozen=True)
class WorkflowSection:
    title: str
    cls: str = "long"
    paths: tuple[tuple[str, tuple[str, ...]], ...] | None = None

    def render(self, graph: ProtocolGraph) -> str:
        return f"<h2>{html.escape(self.title)}</h2>\n" + protocol_graph_panels(
            graph, cls=self.cls, paths=self.paths)


def _sequencing_plan(graph: ProtocolGraph, terminal: str):
    if terminal not in graph.sequencing:
        raise ValueError(f"page readout has no sequencing plan for {terminal!r}")
    return graph.sequencing[terminal]


@dataclass(frozen=True)
class PrimerBindingSection:
    """Display primer landing and computed reads from an exact graph plan."""

    terminal: str
    title: str
    heading: str
    intro_html: str = ""
    after_html: str = ""
    show_read_content: bool = True

    def render(self, graph: ProtocolGraph) -> str:
        plan = _sequencing_plan(graph, self.terminal)
        if plan.library is None or not plan.primers:
            raise ValueError("primer binding view needs exact graph sequencing geometry")
        return (f'<h2 id="seq-primers">{html.escape(self.title)}</h2>\n'
                + (info(self.intro_html) if self.intro_html else "")
                + sp.section(plan.library, plan.primers, heading=self.heading,
                             required_roles=plan.required_roles,
                             read_lengths=plan.read_lengths, read_uses=plan.read_uses,
                             show_read_content=self.show_read_content,
                             unplaced_reads=plan.unplaced_reads,
                             acquisition_unspecified=plan.acquisition_unspecified)
                + sequencing_limitation(plan)
                + self.after_html)


def sequencing_limitation(plan, label: str = "") -> str:
    """Render a sequencing plan's declared limitation, or nothing when it has none."""
    text = (getattr(plan, "limitation", "") or "").strip()
    if not text:
        return ""
    title = "Sequencing limitation" + (f" — {html.escape(label)}" if label else "")
    return info(f"<strong>{title}.</strong> {html.escape(text)}")


@dataclass(frozen=True)
class PlanLimitationSection:
    """Show the limitation recorded on one terminal's graph sequencing plan."""

    terminal: str
    label: str = ""

    def render(self, graph: ProtocolGraph) -> str:
        plan = _sequencing_plan(graph, self.terminal)
        if not plan.limitation.strip():
            raise ValueError(f"sequencing plan for {self.terminal!r} records no limitation")
        return sequencing_limitation(plan, self.label)


@dataclass(frozen=True)
class SourceRead:
    """One source-declared read; display cycles may include dark-cycle recipes."""

    role: str
    cycles: str
    content: str
    spec_role: str
    use: str

    def spec(self):
        fixed = int(self.cycles) if self.cycles.isdigit() else None
        return sp.read_spec(self.spec_role, fixed, self.content.replace(" + ", "; "),
                            self.use)


@dataclass(frozen=True)
class ReadLayoutSection:
    """Show exact primer sites beside a source-declared run recipe."""

    terminal: str
    reads: tuple[SourceRead, ...]
    note: str
    after_html: str = ""

    def render(self, graph: ProtocolGraph) -> str:
        plan = _sequencing_plan(graph, self.terminal)
        if plan.library is None or not plan.primers:
            raise ValueError("read layout view needs exact graph sequencing geometry")
        if sp.verify(plan.library, plan.primers,
                     required_roles=plan.acquired_roles):
            raise ValueError("read layout does not match final library")
        if tuple(read.spec_role for read in self.reads) != plan.acquired_roles:
            raise ValueError("source read rows disagree with graph sequencing roles")
        rows = [(read.role, read.cycles, read.content) for read in self.reads]
        return ('<h2>Read layout</h2>' + sp.diagram(plan.library, plan.primers)
                + table(("Read", "Cycles", "Content"), rows)
                + sp.read_spec_table([read.spec() for read in self.reads], self.note)
                + self.after_html)


@dataclass(frozen=True)
class OffsetPrimerSection:
    """Display a repeated-template primer series without claiming unique read sites."""

    terminal: str
    primers: tuple[tuple[str, str], ...]
    intro: str
    read_content: str
    note: str

    def render(self, graph: ProtocolGraph) -> str:
        plan = _sequencing_plan(graph, self.terminal)
        if not plan.limitation or plan.acquired_roles != ("Sequencing output",):
            raise ValueError("offset-primer view needs a limited sequencing-output plan")
        if not self.primers or any(not name or not sequence for name, sequence in self.primers):
            raise ValueError("offset-primer view needs named primer sequences")
        first = self.primers[0][1]
        if any(sequence != first[offset:] for offset, (_, sequence) in
               enumerate(self.primers)):
            raise ValueError("offset-primer sequences must be successive 5-prime truncations")
        rows = [(name, f"5′-p {sequence}-3′", len(sequence))
                for name, sequence in self.primers]
        return (info(self.intro) + "\n"
                + table(("Primer", "Sequence", "Length"), rows, scroll=False) + "\n"
                + sp.read_spec_table([
                    sp.read_spec("Sequencing output", None, self.read_content, "include")],
                    self.note) + "\n")


def _detail_panel(state, caption: str) -> str:
    if isinstance(state, Construct):
        unknown = tuple(segment.name + "'" for segment in state
                        if segment.is_role_token())
        rows = Scene.duplex(list(state), label=state.name,
                            unpaired=unknown).rows()
    elif isinstance(state, MolecularState):
        rows = state.rows
    elif callable(getattr(state, "rows", None)):
        rows = state.rows()
    elif isinstance(state, (tuple, list)) and all(isinstance(row, Row) for row in state):
        rows = state
    else:
        raise TypeError(f"unsupported graph detail state {type(state).__name__}")
    return panel(rows, cls="long", caption=caption)


def render_annotations(graph: ProtocolGraph, slot: AnnotationSlot) -> str:
    """Render each annotation once in declaration order at its declared page slot."""
    if not isinstance(slot, AnnotationSlot):
        raise ValueError("unknown annotation slot")
    rendered = []
    for annotation in graph.annotations:
        if annotation.slot is not slot:
            continue
        target = html.escape(annotation.target, quote=True)
        title = html.escape(annotation.title)
        if isinstance(annotation, Highlight):
            body = info(f"<strong>{title}.</strong> {html.escape(annotation.explanation)}")
        elif isinstance(annotation, SourceDiscrepancy):
            claims = " ".join(
                f"<li><strong>{html.escape(claim.location)}:</strong> "
                f"{html.escape(claim.claim)}</li>" for claim in annotation.claims)
            body = caveat(
                f"<strong>{title}.</strong><ul>{claims}</ul>"
                f"<p><strong>Modelled:</strong> {html.escape(annotation.model_choice)} "
                f"<strong>Evidence:</strong> {html.escape(annotation.evidence)} "
                f"<strong>Remaining uncertainty:</strong> "
                f"{html.escape(annotation.uncertainty)}</p>")
        else:
            assert isinstance(annotation, DetailView)
            state = annotation.state
            if state is None:
                state = graph.materials[annotation.target].state
            if annotation.style is DetailStyle.FOLDED:
                body = state.folded_panel(annotation.caption)
            else:
                body = _detail_panel(state, annotation.caption)
            body = f"<h2>{title}</h2>\n{body}"
        rendered.append(f'<section data-graph-target="{target}">{body}</section>')
    return "\n".join(rendered)


def assemble_page(graph: ProtocolGraph, parts: Iterable[str | AnnotationSlot |
                                                       OligoSection | WorkflowSection |
                                                       PrimerBindingSection | ReadLayoutSection |
                                                       OffsetPrimerSection |
                                                       PlanLimitationSection]) -> str:
    """Place declared view slots and reject annotations a page forgot to render."""
    output = []
    seen = set()
    for part in parts:
        if isinstance(part, AnnotationSlot):
            if part in seen:
                raise ValueError(f"annotation slot {part.value!r} appears twice")
            seen.add(part)
            output.append(render_annotations(graph, part))
        elif isinstance(part, OligoSection):
            output.append(part.render())
        elif isinstance(part, (WorkflowSection, PrimerBindingSection,
                               ReadLayoutSection, OffsetPrimerSection,
                               PlanLimitationSection)):
            output.append(part.render(graph))
        elif isinstance(part, str):
            output.append(part)
        else:
            raise TypeError("page part must be HTML or a declared view")
    missing = {annotation.slot for annotation in graph.annotations} - seen
    if missing:
        raise ValueError("page omits graph annotation slots: " +
                         ", ".join(sorted(slot.value for slot in missing)))
    return "\n".join(output)
