"""Source-backed specifications for the named-protocol expansion batch.

This module is deliberately data-oriented.  Protocol directories remain separate public
units, while shared constructors make the recurring molecular operations consistent.
"""
from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace
from typing import Callable

import seqprimers as sp
from batch_ngs import nextera_library, truseq_library
from chemdraw import Construct, MolecularState, Scene, Segment, Workflow, feature
from interaction_capture import ReversibleDuplexCrosslink
from rna_library_operations import molecule


def role(name: str, width: int = 18, tag: str | None = None, *, length_bp=None,
         inferred: bool = False, molecular_feature=None) -> Segment:
    token = f"[{length_bp}-bp {name}]" if length_bp else f"[{name}]"
    return Segment(name, token, tag, placeholder=True, length_bp=length_bp,
                   inferred=inferred, feature=molecular_feature)


def two_branch(label: str, a: str, b: str, note: str = "") -> Scene:
    sc = Scene(); sc.strand("branch A", [role(a)], label=label)
    sc.strand("branch B", [role(b)], label=label)
    if note:
        sc.note("branch A", note)
    return sc


def workflow(initial: Scene, steps: tuple[tuple[str, Scene, str], ...]) -> Workflow:
    w = Workflow(MolecularState("Starting material", tuple(initial.rows())))
    for action, state, note in steps:
        w.react(action, state.rows(), name=action, note=note)
    return w


def psoralen_workflow(reagent: str, selection: str) -> Workflow:
    """Crosslink -> fragment and enrich duplexes -> proximity-ligate -> photoreverse.

    This fixed order (duplex enrichment before proximity ligation) is the PARIS and
    SPLASH route.  It is not a general psoralen recipe: a protocol that ligates
    before selecting (e.g. LIGR-seq: circRNA ligase, then RNase R) or installs
    library ends differently must declare its own ordered reactions.
    """
    x = ReversibleDuplexCrosslink(reagent, selection)
    return workflow(molecule("cellular RNA", ("RNA molecules",)), (
        (f"Cross-link paired RNA with {reagent} and long-wave UV", x.duplex(), selection),
        ("Fragment RNA and enrich cross-linked duplexes", x.duplex(), selection),
        ("Proximity-ligate the two RNA partners", x.chimera(),
         f"The covalent ligation is made while {reagent} holds the duplex together."),
        (f"Photoreverse {reagent}, then install sequencing adapters",
         molecule("cDNA library precursor", ("RNA A-derived cDNA", "RNA B-derived cDNA"), joined=True),
         "The chimeric junction, rather than the crosslink, remains in the sequenced molecule."),
    ))

@dataclass(frozen=True)
class Definition:
    """Page text plus a ``library`` shell selector for :func:`get`.

    ``"truseq"`` builds the canonical dual-index TruSeq adapter-ligation endpoint
    with a dA/T ligation junction; ``"nextera"`` the canonical Tn5 endpoint.  These
    are platform defaults, not source facts: do not select them for random-primed,
    template-switching or PCR-added libraries (pass ``final_override`` instead), and
    pass ``inferred_adapters=True`` to :func:`get` when the source does not print
    the adapters.
    """
    title: str
    doi: str
    source_label: str
    summary: str
    library: str
    payload: tuple[str, ...]
    sequencing: str
    modality: str
    caveat: str = ""


DEFS = {
}


def _payload(d: Definition) -> list[Segment]:
    """Role placeholders for ``d.payload``; identifier meaning is read from the names.

    The name matching is a legacy convenience: a caller with source-declared
    features should build its own segments and pass ``final_override``.
    """
    out = []
    for name in d.payload:
        kwargs = {}
        tag = None
        low = name.lower()
        if "umi" in low:
            tag = "umi"; kwargs["feature"] = feature("umi", "umi", "random")
        elif "cell barcode" in low:
            tag = "cbc"; kwargs["feature"] = feature("cell_barcode", "cell_barcode", "unknown")
        elif "sample barcode" in low:
            tag = "cbc"; kwargs["feature"] = feature("sample_index", "sample_index", "unknown")
        elif "spatial hdmi" in low:
            tag = "cbc"; kwargs["feature"] = feature("spatial_hdmi", "spatial_barcode", "random")
        elif "linked-read barcode" in low:
            tag = "cbc"; kwargs["feature"] = feature("linked_read_barcode", "linked_read_barcode", "whitelist", whitelist="stLFR bead barcode set")
        out.append(Segment(name, "X" * (10 if tag else 28), tag, placeholder=True, **kwargs))
    return out


def get(key: str, *, definition: Definition | None = None,
        psoralen: tuple[str, str, str, str] | None = None,
        workflow_fn: Callable[[], Workflow] | None = None,
        final_override=None, primers_override=None, oligos_fn=None,
        sequencing_unavailable: str | None = None,
        sequencing_mechanism: str | None = None,
        motor_adapter_name: str | None = None,
        inferred_adapters: bool = False):
    """Assemble a page from caller-owned source choices and shared library shapes.

    ``inferred_adapters`` marks every platform-shell segment (adapters, index
    windows, junction) as inferred, for sources that do not print them.
    """
    d = definition if definition is not None else DEFS[key]
    primers = ()
    if final_override is not None:
        if isinstance(final_override, tuple):
            final, primers = final_override
        else:
            final = final_override
    elif d.library == "truseq":
        final, primers = truseq_library(_payload(d), d.title + " final library",
                                        inferred_adapters=inferred_adapters)
    elif d.library == "nextera":
        final, primers = nextera_library(_payload(d), d.title + " final library",
                                         inferred_adapters=inferred_adapters)
    elif d.library == "nanopore":
        final = Construct([role(name) for name in d.payload], name=d.title + " final molecule")
    elif d.library == "dnb":
        final = Construct([role("DNBSEQ adapter", inferred=inferred_adapters), *_payload(d),
                           role("DNBSEQ adapter reverse complement", inferred=inferred_adapters)],
                          name=d.title + " final library")
    else:
        final = Construct([role("Illumina-compatible left adapter", inferred=inferred_adapters),
                           *_payload(d),
                           role("Illumina-compatible right adapter", inferred=inferred_adapters)],
                          name=d.title + " final library")
    if primers_override == "standard":
        primers = sp.standard_primers(final, ("Read 1", "Read 2"))
    elif primers_override is not None:
        primers = primers_override
    spec = SimpleNamespace(
        TITLE=d.title, NOTES="01_" + key.replace("-", "_") + ".html",
        SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.source_label}</a>.',
        SUMMARY=d.summary, CAVEAT=d.caveat, FINAL_LIBRARY=final,
        SEQ_PRIMERS=primers, FINAL_CAPTION=d.sequencing,
        SEQUENCING_INTRO=d.sequencing, READ_LENGTHS=None,
        SEQUENCING_ENDING=d.sequencing,
        SEQUENCING_UNAVAILABLE=(sequencing_unavailable if sequencing_unavailable is not None else
                                "The defining source does not expose an exact sequencing-primer sequence; the molecular read entry is shown by role."
                                if not primers else ""),
        workflow=workflow_fn or ((lambda: psoralen_workflow(*psoralen[:2])) if psoralen else None),
        oligos=oligos_fn or (lambda: ()),
    )
    if sequencing_mechanism is not None:
        spec.SEQUENCING_MECHANISM = sequencing_mechanism
    if motor_adapter_name is not None:
        spec.MOTOR_ADAPTER_NAME = motor_adapter_name
    return spec
