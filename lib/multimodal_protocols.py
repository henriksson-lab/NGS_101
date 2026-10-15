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

def _finish_pcr(graph, prefix, template, final, endpoint, reviewed, *, label):
    forward = graph.input(prefix + "_forward", label + " forward PCR primer",
                          MaterialKind.OLIGO)
    reverse = graph.input(prefix + "_reverse", label + " reverse PCR primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input(prefix + "_polymerase", label + " PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr", "Amplify and index the " + label + " library",
                StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="The source defines a primer family or kit-level endpoint rather than one universal amplicon."))
    graph.sequence(library, SequencingPlan(limitation=endpoint))
    reviewed.extend((prefix + "_forward", prefix + "_reverse",
                     prefix + "_polymerase", prefix + "_pcr",
                     prefix + "_library"))
    return library

def _finish_rt_pcr(graph, prefix, template, final, endpoint, reviewed, *, label):
    primer = graph.input(prefix + "_rt_primer", label + " reverse-transcription primer",
                         MaterialKind.OLIGO)
    enzyme = graph.input(prefix + "_rt_enzyme", label + " reverse transcriptase",
                         MaterialKind.CATALYST)
    cdna = graph.product(prefix + "_cdna", label + " cDNA", MaterialKind.MOLECULAR,
                         state=molecule(label + " cDNA", (role("cDNA insert"),)))
    graph.react(prefix + "_rt", "Reverse-transcribe the " + label + " branch",
                StepKind.REVERSE_TRANSCRIPTION,
                ((template, InputRole.TEMPLATE), (primer, InputRole.STRUCTURAL),
                 (enzyme, InputRole.CATALYST)), (cdna,))
    reviewed.extend((prefix + "_rt_primer", prefix + "_rt_enzyme", prefix + "_rt"))
    return _finish_pcr(graph, prefix, cdna, final, endpoint, reviewed, label=label)
