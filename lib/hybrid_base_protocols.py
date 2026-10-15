"""Bounded shared scene and page-spec helpers for conversion and R-loop assays."""
from __future__ import annotations
from dataclasses import dataclass
from types import SimpleNamespace
from chemdraw import Construct, Scene, Segment

def _part(name: str, symbol: str = "X", n: int = 24, **kw) -> Segment:
    return Segment(name, symbol * n, placeholder=True, **kw)


@dataclass(frozen=True)
class RLoopCapture:
    """A three-stranded R-loop plus a capture route and sequenced component."""
    probe: str
    sequenced_component: str
    in_vivo: bool = False

    def __post_init__(self) -> None:
        if self.sequenced_component not in {"DNA", "RNA"}:
            raise ValueError("an R-loop capture must sequence either its DNA or RNA component")
        if not self.probe.strip():
            raise ValueError("an R-loop capture requires an affinity probe")

    def rloop(self) -> Scene:
        sc = Scene()
        sc.strand("template DNA", [_part("DNA in RNA:DNA hybrid")], label="R-loop")
        sc.anneal("hybrid RNA", [_part("RNA in RNA:DNA hybrid", "x")],
                  to="template DNA", pair=("RNA in RNA:DNA hybrid", "DNA in RNA:DNA hybrid"),
                  label="R-loop")
        sc.strand("displaced DNA", [_part("displaced single-stranded DNA", "Y")],
                  label="displaced DNA")
        sc.note("hybrid RNA", f"captured by {self.probe}")
        return sc

    def captured_component(self) -> Scene:
        sc = Scene()
        sc.strand("selected molecule", [_part(f"R-loop {self.sequenced_component} component")],
                  label=f"selected {self.sequenced_component}")
        sc.note("selected molecule", f"{self.probe}; " +
                ("captured in fixed cells" if self.in_vivo else "captured after gentle extraction"))
        return sc


def _role_library(name: str, payload: tuple[str, ...]) -> Construct:
    return Construct([_part("Illumina-compatible left adapter", "L", 18),
                      *[_part(x, "X" if i % 2 == 0 else "Y", 26)
                        for i, x in enumerate(payload)],
                      _part("Illumina-compatible right adapter", "R", 18)], name=name)


def _conversion_scene(states: tuple[tuple[str, str], ...], label: str) -> Scene:
    sc = Scene()
    for i, (origin, state) in enumerate(states):
        sc.strand(origin, [Segment(state, state, "me" if "m" in state else None,
                                   placeholder=True, length_bp=1)], label=f"{origin} → {state}", at=0)
    return sc


def _role_library_scene(name: str, payload: tuple[str, ...], *, ligated: bool = True) -> Scene:
    lib = _role_library(name, payload)
    sc = Scene.duplex(list(lib), label=name)
    if ligated:
        sc.junction("top", "Illumina-compatible left adapter", payload[0], "adapter ligation")
        sc.junction("top", payload[-1], "Illumina-compatible right adapter", "adapter ligation")
    return sc


@dataclass(frozen=True)
class Definition:
    title: str
    doi: str
    source: str
    summary: str
    payload: tuple[str, ...]
    sequencing: str
    notes: str
    truseq: bool = False
    caveat: str = ""


def make_spec(d: Definition, final: Construct, primers: tuple, *,
              oligos=lambda: (), source_detail: str = ""):
    """Assemble common page fields from explicit protocol-owned choices."""
    return SimpleNamespace(
        TITLE=d.title, NOTES=d.notes,
        SOURCE=(f'Defining source: <a href="https://doi.org/{d.doi}">{d.source}</a>.'
                + source_detail),
        SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final,
        SEQ_PRIMERS=primers, FINAL_CAPTION=d.sequencing,
        SEQUENCING_INTRO=d.sequencing, READ_LENGTHS=None,
        SEQUENCING_ENDING=d.sequencing,
        SEQUENCING_UNAVAILABLE=(
            "The authoritative source does not expose the exact sequencing-primer sequence for this library implementation."
            if not primers else ""),
        oligos=oligos)
