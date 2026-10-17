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
                endpoint: str, reviewed: list[str], *, label: str = "library",
                step_label: str | None = None, forward_label: str | None = None,
                reverse_label: str | None = None, polymerase_label: str | None = None,
                forward_state=None, reverse_state=None, primer_source: str | None = None,
                geometry=None,
                limitation: str = "The source defines read roles, but the complete primer geometry is panel-, kit- or implementation-dependent."):
    """Endpoint PCR and sequencing terminal.

    Defaults name roles only and claim neither indexing nor a particular primer;
    callers pass source-named labels, exact printed primer states and geometry.
    """
    forward = graph.input(prefix + "_forward", forward_label or label + " PCR forward primer",
                          MaterialKind.OLIGO, state=forward_state, source=primer_source)
    reverse = graph.input(prefix + "_reverse", reverse_label or label + " PCR reverse primer",
                          MaterialKind.OLIGO, state=reverse_state, source=primer_source)
    polymerase = graph.input(prefix + "_polymerase", polymerase_label or label + " PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr", step_label or "PCR-amplify the " + label, StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=geometry or PCRGeometry.unavailable(
                    final, reason=limitation))
    graph.sequence(library, SequencingPlan(limitation=endpoint))
    reviewed.extend((prefix + "_forward", prefix + "_reverse", prefix + "_polymerase",
                     prefix + "_pcr", prefix + "_library"))
    return library


def _rt_then_pcr(graph: ProtocolGraph, prefix: str, template, final,
                 endpoint: str, reviewed: list[str], *, rt_primer_label: str,
                 label: str = "library", rt_enzyme_label: str = "reverse transcriptase",
                 rt_primer_state=None, between_rt_and_pcr: str | None = None, **pcr):
    """Reverse transcription followed by :func:`_finish_pcr`.

    First-strand cDNA is the PCR template only when ``between_rt_and_pcr`` is
    ``None``.  A caller whose source has further, unmodelled steps (template
    switching, cDNA amplification, IVT, tagmentation, ...) describes them in that
    argument, which becomes an explicit unmodelled boundary step.  Remaining
    keyword arguments go to :func:`_finish_pcr`.
    """
    primer = graph.input(prefix + "_rt_primer", rt_primer_label, MaterialKind.OLIGO,
                         state=rt_primer_state)
    rt = graph.input(prefix + "_rt_enzyme", rt_enzyme_label, MaterialKind.CATALYST)
    cdna = graph.product(prefix + "_cdna", label + " first-strand cDNA",
                         MaterialKind.MOLECULAR,
                         state=molecule(label + " first-strand cDNA", (role("cDNA"),)))
    graph.react(prefix + "_rt", "Reverse-transcribe the " + label + " RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((template, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    reviewed.extend((prefix + "_rt_primer", prefix + "_rt_enzyme", prefix + "_rt"))
    pcr_template = cdna
    if between_rt_and_pcr is not None:
        pcr_template = graph.product(
            prefix + "_pcr_precursor", label + " PCR template (steps not modelled)",
            MaterialKind.MOLECULAR,
            state=molecule(label + " PCR template", (role("library precursor"),)))
        graph.react(prefix + "_unmodelled",
                    "Source steps between reverse transcription and endpoint PCR (not modelled)",
                    StepKind.MAGIC, ((cdna, InputRole.TEMPLATE),), (pcr_template,),
                    reason=between_rt_and_pcr)
        reviewed.append(prefix + "_unmodelled")
    return _finish_pcr(graph, prefix, pcr_template, final, endpoint, reviewed,
                       label=label, **pcr)
