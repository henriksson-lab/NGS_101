"""Structured evidence for the reagent-end geometry of adapter ligations.

The sealed library cannot by itself distinguish a blunt, cohesive, splinted or simply
undisclosed ligation substrate.  Protocol constructors therefore pass one of these
objects to :meth:`chemdraw.Workflow.ligate_adapter` before producing the sealed state.
"""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import Row, Scene, Segment


@dataclass(frozen=True)
class AdapterLigationGeometry:
    availability: str
    kind: str
    scenes: tuple[Scene, ...]
    statement: str

    def __post_init__(self) -> None:
        if self.availability not in ("known", "unavailable"):
            raise ValueError("adapter geometry must be known or unavailable")
        if self.kind not in ("cohesive", "blunt", "splinted", "single-strand",
                             "unavailable"):
            raise ValueError("unknown adapter-ligation geometry kind")
        if not self.statement.strip():
            raise ValueError("adapter geometry needs a source-scoped statement")
        adapter_ends = [e for sc in self.scenes for e in sc.ligation_ends
                        if e.owner == "adapter"]
        if self.availability == "known" and not any(e.kind == self.kind for e in adapter_ends):
            raise ValueError("known adapter geometry must declare its adapter end")
        if self.availability == "unavailable" and adapter_ends:
            raise ValueError("unavailable adapter geometry cannot declare a known end")
        if (self.availability == "unavailable") != (self.kind == "unavailable"):
            raise ValueError("only unavailable evidence may use unavailable geometry")

    @classmethod
    def known(cls, kind: str, *scenes: Scene,
              statement: str) -> "AdapterLigationGeometry":
        if not scenes:
            raise ValueError("known adapter geometry needs a molecular scene")
        return cls("known", kind, tuple(scenes), statement)

    @classmethod
    def unavailable(cls, reason: str) -> "AdapterLigationGeometry":
        sc = Scene()
        sc.strand("adapter end", [Segment("adapter overhang not reported",
                                          "[adapter overhang not reported]",
                                          placeholder=True)],
                  label="source limitation")
        sc.note("adapter end", reason)
        return cls("unavailable", "unavailable", (sc,), reason)

    def rows(self):
        rows = []
        for i, sc in enumerate(self.scenes):
            if i:
                rows.append(Row())
            rows.extend(sc.rows())
        return rows


def reported_adapter_end(name: str, *, kind: str, statement: str,
                         polarity: str | None = None,
                         terminal_sequence: str | None = None,
                         body_sequence: str | None = None,
                         substrate_scenes: tuple[Scene, ...] = ()) -> AdapterLigationGeometry:
    """Construct a source-known terminal adapter end without inventing its body.

    This is for sources that disclose the ligating end (for example a 3′ dT) but keep
    the rest of a commercial adapter proprietary.  Exact full duplexes should instead be
    built and paired explicitly with :class:`Scene`.
    """
    if kind in ("cohesive", "splinted", "single-strand") and not terminal_sequence:
        raise ValueError(f"{kind} adapter end needs its reported terminal sequence")
    body = Segment(f"{name} body",
                   body_sequence or "[adapter body not disclosed]",
                   placeholder=body_sequence is None)
    end = Segment(f"{name} ligating end", terminal_sequence or "[blunt end]",
                  placeholder=terminal_sequence is None)
    parts = ([end, body] if polarity == "5-prime" else [body, end])
    sc = Scene(); sc.strand("adapter", parts, label=name)
    sc.ligation_end("adapter", end.name, kind=kind, polarity=polarity,
                    owner="adapter")
    return AdapterLigationGeometry.known(
        kind, sc, *substrate_scenes, statement=statement)
