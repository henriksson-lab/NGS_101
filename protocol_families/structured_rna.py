"""Bounded operations for structured-RNA library preparation."""
from __future__ import annotations
from dataclasses import dataclass
from chemdraw import Construct, Scene, Segment
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role

def exact(name: str, seq: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name, seq, tag, **kw)


@dataclass(frozen=True)
class StructuredRNAInput:
    name: str
    end_state: str
    demethylation: str = "none"
    preserve_modification_signatures: bool = True

    def __post_init__(self) -> None:
        if self.demethylation != "none" and self.preserve_modification_signatures:
            raise ValueError("demethylation cannot preserve the removed modification signatures")

    def scene(self) -> Scene:
        sc = molecule(self.name, (role("structured RNA body"), role(self.end_state)))
        sc.note("molecule", f"demethylation: {self.demethylation}")
        return sc


def _finish_rt_pcr(graph: ProtocolGraph, prefix: str, template, final: Construct,
                   endpoint: str, reviewed: list[str], *, rt_primer_label: str):
    rt_primer = graph.input(prefix + "_rt_primer", rt_primer_label, MaterialKind.OLIGO)
    reverse_transcriptase = graph.input(prefix + "_reverse_transcriptase", "reverse transcriptase",
                                       MaterialKind.CATALYST)
    cdna = graph.product(prefix + "_cdna", "structured-RNA-derived cDNA",
                         MaterialKind.MOLECULAR,
                         state=molecule("structured-RNA-derived cDNA", (
                             role("structured RNA cDNA"),)))
    graph.react(prefix + "_rt", "Reverse-transcribe the adapted RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((template, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (cdna,))
    reviewed.extend((prefix + "_rt_primer", prefix + "_reverse_transcriptase",
                     prefix + "_rt"))
    return _finish_pcr(graph, prefix, cdna, final, endpoint, reviewed)


def _finish_pcr(graph: ProtocolGraph, prefix: str, template, final: Construct,
                endpoint: str, reviewed: list[str], *, sequence: bool = True):
    forward = graph.input(prefix + "_pcr_forward", "library PCR forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input(prefix + "_pcr_reverse", "indexed library PCR reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input(prefix + "_pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", "indexed structured-RNA library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr", "Index and amplify the library", StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final,
                    reason="The source defines the PCR oligo roles but not a distinct sequencing-primer landing geometry on this role-level construct."))
    if sequence:
        graph.sequence(library, SequencingPlan(limitation=endpoint))
    reviewed.extend((
        prefix + "_pcr_forward", prefix + "_pcr_reverse", prefix + "_pcr_polymerase",
        prefix + "_pcr", prefix + "_library",
    ))
    return library
