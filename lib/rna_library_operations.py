"""Legacy RNA-library recipes awaiting source-specific caller migration.

These helpers choose reaction order, unnamed reagents, role lengths, review scope
and source-absence wording. Supplying a final construct does not make those choices
neutral or source-verified. Do not add new callers; declare the source-specific
route locally and use shared molecular operations only with explicit inputs.
See catalogue/chemistry_decoupling.tsv.
"""
from __future__ import annotations

from adapter_ligation import AdapterLigationGeometry
from chemdraw import Construct, Scene, Segment
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)
from protocol_model import molecule as role_molecule, role

def molecule(label: str, names: tuple[str, ...], *, joined: bool = False,
             note: str = "", junctions: tuple[tuple[str, str, str], ...] = ()) -> Scene:
    parts = [Segment(name, chr(88 + (i % 2)) * 18, placeholder=True)
             for i, name in enumerate(names)]
    sc = Scene(); sc.strand("molecule", parts, label=label)
    # ``joined`` means only that the labels describe one covalent molecule.  A barcode
    # copied from a primer is just as joined as a ligated adapter, but it is not a
    # ligation site.  Only callers naming an exact boundary may draw the ** marker.
    for left, right, label_text in junctions:
        sc.junction("molecule", left, right, label_text)
    if note:
        sc.note("molecule", note)
    sc.labels("molecule")
    return sc

def finish_pcr(graph, prefix, template, final, primers, endpoint, reviewed, *, label,
               step_label=None, forward_label=None, reverse_label=None,
               polymerase_label=None, forward_state=None, reverse_state=None,
               primer_source=None, primers_inferred=False, geometry=None,
               unavailable_reason=None):
    """Legacy endpoint PCR.  Default labels name roles only, never an index claim.

    Callers name source primers (e.g. an indexed reverse primer), pass printed
    primer states and an exact ``geometry`` or a source-specific
    ``unavailable_reason``; ``primers_inferred`` marks unprinted primer inputs.
    """
    forward = graph.input(prefix + "_forward", forward_label or label + " forward PCR primer",
                          MaterialKind.OLIGO, state=forward_state, source=primer_source,
                          inferred=primers_inferred)
    reverse = graph.input(prefix + "_reverse", reverse_label or label + " reverse PCR primer",
                          MaterialKind.OLIGO, state=reverse_state, source=primer_source,
                          inferred=primers_inferred)
    polymerase = graph.input(prefix + "_polymerase", polymerase_label or label + " PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", label + " sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr", step_label or "Amplify the " + label + " library",
                StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=geometry or PCRGeometry.unavailable(
                    final, reason=unavailable_reason or "The defining source establishes the library roles but not one universal source-exact amplicon geometry."))
    plan = (SequencingPlan(library=final, primers=tuple(primers)) if primers
            else SequencingPlan(limitation=endpoint))
    graph.sequence(library, plan)
    reviewed.extend((prefix + "_forward", prefix + "_reverse", prefix + "_polymerase",
                     prefix + "_pcr", prefix + "_library"))
    return library


def finish_rt_pcr(graph, prefix, template, final, primers, endpoint, reviewed, *, label,
                  **pcr_options):
    rt_primer = graph.input(prefix + "_rt_primer", label + " reverse-transcription primer",
                            MaterialKind.OLIGO)
    rt = graph.input(prefix + "_rt_enzyme", label + " reverse transcriptase",
                     MaterialKind.CATALYST)
    cdna = graph.product(prefix + "_cdna", label + " cDNA", MaterialKind.MOLECULAR,
                         state=molecule(label + " cDNA", ("cDNA insert",)))
    graph.react(prefix + "_rt", "Reverse-transcribe the " + label + " product",
                StepKind.REVERSE_TRANSCRIPTION,
                ((template, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (cdna,))
    reviewed.extend((prefix + "_rt_primer", prefix + "_rt_enzyme", prefix + "_rt"))
    return finish_pcr(graph, prefix, cdna, final, primers, endpoint, reviewed, label=label,
                      **pcr_options)

def ligate_rna_adapter(graph: ProtocolGraph, prefix: str, template, label: str,
                        reviewed: list[str]):
    adapter = graph.input(prefix + "_adapter", label + " RNA adapter", MaterialKind.ADAPTER)
    ligase = graph.input(prefix + "_ligase", label + " RNA ligase", MaterialKind.CATALYST)
    product = graph.product(prefix + "_product", label + " adapter-ligated RNA",
                            MaterialKind.MOLECULAR,
                            state=role_molecule(label + " adapter-ligated RNA", (
                                role(label + " adapter"), role("selected RNA"))))
    graph.react(prefix + "_ligation", "Ligate the " + label + " RNA adapter",
                StepKind.LIGATION,
                ((template, InputRole.TEMPLATE), (adapter, InputRole.STRUCTURAL),
                 (ligase, InputRole.CATALYST)), (product,),
                geometry=AdapterLigationGeometry.unavailable(
                    "The defining source establishes this RNA-adapter role but not one immutable ligating-end sequence for the page model."))
    reviewed.extend((prefix + "_adapter", prefix + "_ligase", prefix + "_ligation"))
    return product


def finish_cdna(graph: ProtocolGraph, prefix: str, template, final: Construct,
                 endpoint: str, reviewed: list[str], *, rt_label: str = "adapter RT primer",
                 indexed: bool = True):
    rt_primer = graph.input(prefix + "_rt_primer", rt_label, MaterialKind.OLIGO)
    reverse_transcriptase = graph.input(prefix + "_reverse_transcriptase", "reverse transcriptase",
                                       MaterialKind.CATALYST)
    cdna = graph.product(prefix + "_cdna", "assay-derived first-strand cDNA",
                         MaterialKind.MOLECULAR,
                         state=role_molecule("assay-derived cDNA", (role("selected RNA cDNA"),)))
    graph.react(prefix + "_rt", "Reverse-transcribe the selected RNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((template, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (cdna,))
    forward = graph.input(prefix + "_pcr_forward", "library PCR forward primer",
                          MaterialKind.OLIGO)
    reverse = graph.input(prefix + "_pcr_reverse", "indexed library PCR reverse primer" if indexed else "library PCR reverse primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input(prefix + "_pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product(prefix + "_library", "indexed sequencing library" if indexed else "sequencing library",
                            MaterialKind.LIBRARY, state=final)
    graph.react(prefix + "_pcr", "Index and amplify the cDNA library" if indexed else "Amplify the cDNA library", StepKind.PCR,
                ((cdna, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final,
                    reason="The source establishes read direction and library roles but not one immutable complete primer geometry."))
    graph.sequence(library, SequencingPlan(limitation=endpoint))
    reviewed.extend((
        prefix + "_rt_primer", prefix + "_reverse_transcriptase", prefix + "_rt",
        prefix + "_pcr_forward", prefix + "_pcr_reverse", prefix + "_pcr_polymerase",
        prefix + "_pcr", prefix + "_library",
    ))
    return library
