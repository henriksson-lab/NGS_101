"""Structured evidence for the reagent-end geometry of adapter ligations.

The sealed library cannot by itself distinguish a blunt, cohesive, splinted or simply
undisclosed ligation substrate.  Protocol constructors therefore pass one of these
objects to :meth:`chemdraw.Workflow.ligate_adapter` before producing the sealed state.
"""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import MolecularFeature, Row, Scene, Segment, revcomp


@dataclass(frozen=True)
class CohesiveDuplex:
    """A duplex reconstructed from two source-written 5′→3′ oligos.

    The common case for split-pool barcode plates is a paired core flanked by a
    5′ overhang on each strand.  Callers provide the oligos exactly as printed;
    this constructor finds the longest exact antiparallel overlap and rejects
    internal, 3′-overhang, blunt or ambiguous arrangements.  The resulting
    :class:`Scene` therefore cannot silently draw an incompatible sticky end.
    """

    name: str
    top_5to3: str
    bottom_5to3: str
    top_overhang: str
    paired_core: str
    bottom_overhang: str
    scene: Scene


def cohesive_duplex(name: str, top_5to3: str, bottom_5to3: str, *,
                     core_name: str = "paired core",
                     core_feature: MolecularFeature | None = None,
                     top_mod5: str = "", bottom_mod5: str = "",
                     bottom_mod3: str = "") -> CohesiveDuplex:
    """Build a checked two-sided 5′-cohesive duplex from written oligos.

    One side may be blunt (the common flanking adapters in an assembly often
    expose only one cohesive end), but at least one 5′ overhang is required.
    """
    top = top_5to3.upper().replace(" ", "")
    bottom = bottom_5to3.upper().replace(" ", "")
    if not top or not bottom or set(top + bottom) - set("ACGT"):
        raise ValueError("cohesive-duplex oligos must be non-empty DNA sequences")
    bottom_drawn = revcomp(bottom)
    candidates = []
    for offset in range(0, len(top) + 1):
        end = min(len(top), offset + len(bottom_drawn))
        if end <= offset:
            continue
        overlap = top[offset:end]
        if overlap == bottom_drawn[:len(overlap)] and end == len(top):
            # A 5′ overhang on the bottom strand appears to the right when drawn.
            candidates.append((len(overlap), offset, bottom_drawn[len(overlap):]))
    if not candidates:
        raise ValueError(f"{name}: oligos do not form a terminal antiparallel duplex")
    longest = max(n for n, _, _ in candidates)
    best = [(n, offset, right) for n, offset, right in candidates if n == longest]
    if len(best) != 1:
        raise ValueError(f"{name}: ambiguous longest antiparallel overlap")
    _, offset, right = best[0]
    left, core = top[:offset], top[offset:]
    if not left and not right:
        raise ValueError(f"{name}: duplex is blunt, not cohesive")

    top_parts = []
    if left:
        top_parts.append(Segment(f"{name} left 5′ overhang", left))
    top_parts.append(Segment(core_name, core, feature=core_feature))
    bottom_parts = []
    if right:
        bottom_parts.append(Segment(f"{name} right 5′ overhang", revcomp(right)))
    bottom_parts.append(Segment(f"{core_name} complement", revcomp(core),
                                feature=core_feature))
    sc = Scene()
    sc.strand("top", top_parts, label=f"{name} top", mod5=top_mod5)
    sc.anneal("bottom", bottom_parts, to="top",
              pair=(f"{core_name} complement", core_name), label=f"{name} bottom",
              mod5=bottom_mod5, mod3=bottom_mod3)
    if left:
        sc.overhang("top", f"{name} left 5′ overhang", polarity="5-prime",
                    text=f"5′-{left} cohesive end")
    if right:
        sc.overhang("bottom", f"{name} right 5′ overhang", polarity="5-prime",
                    text=f"5′-{revcomp(right)} cohesive end")
    sc.labels("top")
    return CohesiveDuplex(name, top, bottom, left, core, revcomp(right), sc)


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
                   body_sequence or "[adapter body not declared in this model]",
                   placeholder=body_sequence is None)
    end = Segment(f"{name} ligating end", terminal_sequence or "[blunt end]",
                  placeholder=terminal_sequence is None)
    parts = ([end, body] if polarity == "5-prime" else [body, end])
    sc = Scene(); sc.strand("adapter", parts, label=name)
    sc.ligation_end("adapter", end.name, kind=kind, polarity=polarity,
                    owner="adapter")
    return AdapterLigationGeometry.known(
        kind, sc, *substrate_scenes, statement=statement)
