"""Bounded role and amplification helpers for extended protocols."""
from __future__ import annotations

from dataclasses import dataclass
from chemdraw import feature
from pcr import PCRGeometry
from protocol_graph import InputRole, MaterialKind, ProtocolGraph, SequencingPlan, StepKind
from protocol_model import molecule, role


def ident(name, identifier, kind, *, length=None, tag="cbc"):
    return role(name, tag=tag, length_bp=length,
                molecular_feature=feature(identifier, kind, "unknown"))


CELL = ident("cell barcode", "cell_barcode", "cell_barcode")
UMI = ident("UMI", "umi", "umi", tag="umi")
SAMPLE = ident("sample index", "sample_index", "sample_index", length=8)
TARGET = ident("target barcode", "target_barcode", "feature_barcode")
GUIDE = ident("guide identity", "guide_barcode", "guide_barcode")
RT_WELL = ident("RT-well barcode", "rt_well", "cell_barcode")
LIG_WELL = ident("ligation-well barcode", "ligation_well", "cell_barcode")
PCR_WELL = ident("PCR-well barcode", "pcr_well", "cell_barcode")


@dataclass(frozen=True)
class Definition:
    title: str; doi: str; citation: str; summary: str; endpoint: str; notes: str


def _finish_pcr(graph: ProtocolGraph, prefix: str, template, final,
                endpoint: str, reviewed: list[str], *, label: str = "library"):
    forward = graph.input(prefix + "_forward", label + " PCR forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input(prefix + "_reverse", label + " indexed PCR reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input(prefix + "_polymerase", label + " PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr", "Index and amplify the " + label, StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final,
                    reason="The source defines read roles, but the complete primer geometry is panel-, kit- or implementation-dependent."))
    graph.sequence(library, SequencingPlan(limitation=endpoint))
    reviewed.extend((prefix + "_forward", prefix + "_reverse", prefix + "_polymerase",
                     prefix + "_pcr", prefix + "_library"))
    return library


def _rt_then_pcr(graph: ProtocolGraph, prefix: str, template, final,
                 endpoint: str, reviewed: list[str], *, rt_primer_label: str,
                 label: str = "library"):
    primer = graph.input(prefix + "_rt_primer", rt_primer_label, MaterialKind.OLIGO)
    rt = graph.input(prefix + "_rt_enzyme", "reverse transcriptase", MaterialKind.CATALYST)
    cdna = graph.product(prefix + "_cdna", label + " first-strand cDNA",
                         MaterialKind.MOLECULAR,
                         state=molecule(label + " first-strand cDNA", (role("cDNA"),)))
    graph.react(prefix + "_rt", "Reverse-transcribe the " + label + " RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((template, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    reviewed.extend((prefix + "_rt_primer", prefix + "_rt_enzyme", prefix + "_rt"))
    return _finish_pcr(graph, prefix, cdna, final, endpoint, reviewed, label=label)
