"""Neutral nascent-RNA capture roles and page-spec assembly."""
from __future__ import annotations
from dataclasses import dataclass
from types import SimpleNamespace
from chemdraw import Construct, Scene, Segment, feature
from protocol_model import molecule, role
from rna_library_operations import finish_cdna, ligate_rna_adapter
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind

@dataclass(frozen=True)
class NascentEndCapture:
    """A source of nascent RNA and the end whose coordinate is informative."""
    source: str
    selected_end: str
    cap_selected: bool = False

    def __post_init__(self) -> None:
        if self.selected_end not in {"5-prime", "3-prime", "both"}:
            raise ValueError("nascent-RNA capture must declare 5-prime, 3-prime or both")

    def RNA(self, name: str = "selected nascent RNA") -> Scene:
        parts = [role("5-prime cap" if self.cap_selected else "native 5-prime end"),
                 role("nascent RNA body"), role("polymerase-proximal 3-prime end")]
        sc = molecule(name, tuple(parts))
        sc.note("molecule", f"{self.source}; informative end: {self.selected_end}")
        return sc


@dataclass(frozen=True)
class FourSUReadout:
    """A validated 4sU route: enrichment and base recoding are independent choices."""
    chemistry: str
    output: str
    enrichment: bool

    def __post_init__(self) -> None:
        if self.output not in {"selected 4sU-RNA", "U-to-C mismatch", "native RNA ends"}:
            raise ValueError("unknown 4sU readout")
        if self.output == "U-to-C mismatch" and not self.chemistry:
            raise ValueError("base recoding requires a named chemistry")

    def labelled(self) -> Scene:
        return molecule("pulse-labelled RNA", (role("pre-existing RNA"),
                        Segment("incorporated 4sU", "4sU", placeholder=True, length_bp=1),
                        role("new RNA")))

    def product(self) -> Scene:
        if self.output == "U-to-C mismatch":
            mark = Segment("recoded 4sU read as C", "C", "me", placeholder=True,
                           length_bp=1)
        else:
            mark = Segment(self.output, "U", "me", placeholder=True, length_bp=1)
        sc = molecule(self.output, (role("RNA before label"), mark,
                                    role("RNA after label")))
        sc.note("molecule", self.chemistry)
        return sc


SAMPLE_INDEX = role("sample index", tag="cbc", length_bp=8,
                    molecular_feature=feature("sample_index", "sample_index", "unknown"))


def role_library(name: str, payload: tuple[Segment, ...], *, indexed: bool = True,
                 platform: str = "Illumina-compatible") -> Construct:
    parts = [role(f"{platform} left end")]
    if indexed:
        parts.append(SAMPLE_INDEX)
    parts.extend(payload)
    parts.append(role(f"{platform} right end"))
    return Construct(parts, name=name)


@dataclass(frozen=True)
class Definition:
    title: str
    doi: str
    citation: str
    summary: str
    endpoint: str
    notes: str
    caveat: str = ""


def make_spec(d: Definition, final: Construct, *, unavailable: str,
              oligos=lambda: ()):
    """Assemble page fields from choices declared by the owning protocol."""
    return SimpleNamespace(
        TITLE=d.title, NOTES=d.notes,
        SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
        SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final, SEQ_PRIMERS=(),
        FINAL_CAPTION=d.endpoint, SEQUENCING_INTRO=d.endpoint, READ_LENGTHS=None,
        SEQUENCING_ENDING=d.endpoint, SEQUENCING_UNAVAILABLE=unavailable,
        oligos=oligos)
