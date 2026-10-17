"""Shared molecular models for single-cell multimodal identifier transfer.

The central invariant is identifier provenance: every feature says whether it identifies
a cell, molecule, sample/pool, antibody or genomic target.  Branching workflows keep the
same cell identity while allowing physically distinct sequencing libraries.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace

from chemdraw import feature
from circular import CircularizationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import chain, molecule, role, role_library

def ident(name, identifier, kind, *, length=None, tag="cbc"):
    return role(name, tag=tag, length_bp=length,
                molecular_feature=feature(identifier, kind, "unknown"))

CELL = ident("cell barcode", "cell_barcode", "cell_barcode")
UMI = ident("UMI", "umi", "umi", tag="umi")
SAMPLE = ident("sample index", "sample_index", "sample_index", length=8)
POOL = ident("staining-pool barcode", "pool_barcode", "sample_index", length=10)
ANTIBODY = ident("antibody barcode", "antibody_barcode", "feature_barcode")
GUIDE = ident("12-nt internal guide barcode", "guide_ibar", "guide_barcode", length=12)

@dataclass(frozen=True)
class Definition:
    title: str; doi: str; citation: str; summary: str; endpoint: str; notes: str

def _finish_pcr(graph, prefix, template, final, endpoint, reviewed, *, label,
                step_label=None, forward_label=None, reverse_label=None,
                polymerase_label=None, forward_state=None, reverse_state=None,
                primer_source=None, geometry=None,
                limitation="The source defines a primer family or kit-level endpoint rather than one universal amplicon."):
    """Endpoint PCR and sequencing terminal for one branch.

    Defaults are deliberately unspecific stubs: they name roles, not source chemistry,
    and do not claim indexing.  Callers pass source-named labels, exact printed primer
    states (``forward_state`` / ``reverse_state``, located by ``primer_source``) and an
    exact ``geometry`` when their source supports them.
    """
    forward = graph.input(prefix + "_forward",
                          forward_label or f"unspecified forward PCR primer ({label})",
                          MaterialKind.OLIGO, state=forward_state, source=primer_source)
    reverse = graph.input(prefix + "_reverse",
                          reverse_label or f"unspecified reverse PCR primer ({label})",
                          MaterialKind.OLIGO, state=reverse_state, source=primer_source)
    polymerase = graph.input(prefix + "_polymerase",
                             polymerase_label or f"unspecified PCR polymerase ({label})",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr", step_label or "PCR-amplify the " + label + " library",
                StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=geometry or PCRGeometry.unavailable(
                    final, reason=limitation))
    graph.sequence(library, SequencingPlan(limitation=endpoint))
    reviewed.extend((prefix + "_forward", prefix + "_reverse",
                     prefix + "_polymerase", prefix + "_pcr",
                     prefix + "_library"))
    return library

def _finish_rt_pcr(graph, prefix, template, final, endpoint, reviewed, *, label,
                   rt_primer_label=None, rt_enzyme_label=None, rt_primer_state=None,
                   rt_step_label=None, between_rt_and_pcr=None, **pcr):
    """Reverse transcription followed by :func:`_finish_pcr`.

    The RT defaults are unspecific stubs.  The first-strand product goes straight into
    the endpoint PCR only when ``between_rt_and_pcr`` is ``None``; a protocol whose
    source has further steps (template switching, cDNA amplification, IVT,
    tagmentation, ...) that it does not model passes their description as that
    argument, which becomes an explicit unmodelled boundary step.  Remaining keyword
    arguments go to :func:`_finish_pcr`.
    """
    primer = graph.input(prefix + "_rt_primer",
                         rt_primer_label or f"unspecified reverse-transcription primer ({label})",
                         MaterialKind.OLIGO, state=rt_primer_state)
    enzyme = graph.input(prefix + "_rt_enzyme",
                         rt_enzyme_label or f"unspecified reverse transcriptase ({label})",
                         MaterialKind.CATALYST)
    cdna = graph.product(prefix + "_cdna", label + " cDNA", MaterialKind.MOLECULAR,
                         state=molecule(label + " cDNA", (role("cDNA insert"),)))
    graph.react(prefix + "_rt", rt_step_label or "Reverse-transcribe the " + label + " branch",
                StepKind.REVERSE_TRANSCRIPTION,
                ((template, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
                 (enzyme, InputRole.CATALYST)), (cdna,))
    reviewed.extend((prefix + "_rt_primer", prefix + "_rt_enzyme", prefix + "_rt"))
    pcr_template = _unmodelled_boundary(graph, prefix, cdna, label, between_rt_and_pcr, reviewed)
    return _finish_pcr(graph, prefix, pcr_template, final, endpoint, reviewed,
                       label=label, **pcr)

def _unmodelled_boundary(graph, prefix, template, label, reason, reviewed):
    """Return ``template`` or, when ``reason`` is given, an explicit unmodelled step's product."""
    if reason is None:
        return template
    precursor = graph.product(prefix + "_pcr_precursor", label + " PCR template (steps not modelled)",
                              MaterialKind.MOLECULAR,
                              state=molecule(label + " PCR template", (role("library precursor"),)))
    graph.react(prefix + "_unmodelled", "Source steps between reverse transcription and endpoint PCR (not modelled)",
                StepKind.MAGIC, ((template, InputRole.TEMPLATE),), (precursor,), reason=reason)
    reviewed.append(prefix + "_unmodelled")
    return precursor
