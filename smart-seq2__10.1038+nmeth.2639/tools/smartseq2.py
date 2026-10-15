"""Canonical SMART-seq2 oligos and reusable cDNA/PCR geometry.

Protocol pages that explicitly inherit SMART-seq2 should construct the branch here
instead of restating its oligos or calling them unavailable.
"""
from __future__ import annotations

from chemdraw import Construct, Scene, Segment, feature, revcomp
from illumina import P5, P7_RC
import nextera as nx
from pcr import PCRGeometry, PCRPrimer
import rt
import seqprimers as sp
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

DT_LEN = 30
DT_ANCHOR = "VN"
DT_LINKER = "AC"
TSO_LINKER = "ACAT"
OLIGO_DT = rt.oligo_dt(DT_LEN, DT_ANCHOR, rt.SMART_HANDLE + DT_LINKER)
TSO = rt.SMART_HANDLE + TSO_LINKER + rt.TSO_G_TAIL_LNA
ISPCR = rt.SMART_HANDLE
SEQ_PRIMERS = tuple(sp.NEXTERA[k] for k in ("R1", "I1", "I2", "R2"))
DOI = "10.1038/nmeth.2639"
CAVEAT = ("<div class='caveat'><b>Two published TSO linkers.</b> Picelli et al. "
          "2013 Supplementary Table 4 names <code>rGrG+G</code> as the Smart-seq2 "
          "template-switch oligo, and Supplementary Table 2 prints that oligo as "
          "the 23-nt handle + <code>AC</code> + <code>rGrG+G</code>. The Nature "
          "Protocols version of the protocol instead carries <code>ACAT</code> "
          "before the G tail &mdash; the form transcribed in Hagemann-Jensen et al. "
          "2020 Supplementary Table 1 (&ldquo;Smart-seq2 as of Nature Protocol&rdquo;) "
          "and the one modelled here. The two bases the longer form adds beyond the "
          "defining paper are drawn dotted as inferred.</div>")
TEMPLATE_SWITCH_CAPTION = (
    "INFERRED — Template-switch at the transcript 5-prime end; the TSO "
    "linker's final AT is the Nature Protocols form, not the defining paper's")


def defining_tso_linker() -> list[Segment]:
    """The 2013 AC bases and later AT extension, with their distinct evidence."""
    if not TSO_LINKER.startswith(DT_LINKER):
        raise ValueError("the later TSO linker must extend the published AC linker")
    return [Segment("AC", DT_LINKER),
            Segment("AT", TSO_LINKER[len(DT_LINKER):], inferred=True)]


def drawn_tso() -> list[Segment]:
    """The 2013 source boundary as segments for the template-switch state."""
    return [Segment("H", rt.SMART_HANDLE, "tso"),
            *defining_tso_linker(), Segment("rGrG", "GGG", "tso")]


def drawn_preamp_template(body_nt: int, poly_a_nt: int) -> Construct:
    """Source-marked cDNA state used by the Smart-seq2 reaction graph."""
    if body_nt < 1 or poly_a_nt != DT_LEN:
        raise ValueError("drawn cDNA needs a body and the published dT30 tail")
    link = defining_tso_linker()
    return Construct([
        Segment("H5", rt.SMART_HANDLE, "tso"),
        Segment("TSO linker", link[0].top),
        Segment("TSO linker AT", link[1].top, inferred=True),
        Segment("GGG", "GGG", "tso"),
        Segment("body", "X" * body_nt, placeholder=True),
        Segment("B", "B", placeholder=True),
        Segment("polyA", "A" * poly_a_nt),
        Segment("GT", revcomp(DT_LINKER)),
        Segment("H3", revcomp(rt.SMART_HANDLE), "tso"),
    ], name="SMART-seq full-length cDNA before ISPCR")


def final_library() -> Construct:
    """Heterologous Nextera fragment after index PCR; no internal UMI."""
    return Construct([
        Segment("Illumina P5", P5, "p5"),
        Segment("i5", "N" * 8, placeholder=True,
                feature=feature("sample_index_i5", "sample_index", "unknown")),
        Segment("s5", nx.S5, "s5"), Segment("ME", nx.ME, "me"),
        Segment("cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
        Segment("ME ", nx.ME_RC, "me"), Segment("s7", nx.S7_RC, "s7"),
        Segment("i7", "N" * 8, placeholder=True,
                feature=feature("sample_index_i7", "sample_index", "unknown")),
        Segment("Illumina P7", P7_RC, "p7"),
    ], name="SMART-seq2 library")


def protocol_graph() -> ProtocolGraph:
    # The published predecessor supplies molecular drawing scenes used by this method.
    import sys
    from pathlib import Path
    family = Path(__file__).resolve().parents[2] / "smart-seq__10.1038+nbt.2282" / "tools"
    sys.path.insert(0, str(family))
    from smartseq_states import BODY, POLYA, gap_filled, mrna, ss2_rt, tagmented
    protocol, doi = "SMART-seq2", DOI
    source = f"Defining source: doi:{doi}"
    graph = ProtocolGraph(protocol, authority=source)
    initial = Scene(); initial.strand("mRNA", mrna())
    mrna = graph.input("mrna", "single-cell poly(A) RNA", MaterialKind.MOLECULAR,
                       state=initial, primary=True)

    priming, first_scene, switch_scene = ss2_rt(protocol)
    rt_label = "anchored oligo-dTVN with SMART handle"
    tso_label = "SMART-seq2 LNA template-switch oligo"
    preamp_template = drawn_preamp_template(BODY, POLYA)

    rt_primer = graph.input("rt_primer", rt_label, MaterialKind.OLIGO,
                            state=priming)
    reverse_transcriptase = graph.input(
        "reverse_transcriptase", "MMLV-derived template-switching reverse transcriptase",
        MaterialKind.CATALYST)
    first = graph.product("first_strand", "first-strand cDNA with untemplated CCC",
                          MaterialKind.MOLECULAR, state=first_scene)
    graph.react("reverse_transcription",
                "Reverse-transcribe from the anchored poly(dT) primer",
                StepKind.REVERSE_TRANSCRIPTION,
                ((mrna, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (first,))

    tso = graph.input("tso", tso_label, MaterialKind.OLIGO, state=switch_scene)
    switched = graph.product("switched", "template-switched full-length cDNA",
                             MaterialKind.MOLECULAR, state=preamp_template)
    switch_label = TEMPLATE_SWITCH_CAPTION
    graph.react("template_switch", switch_label,
                StepKind.TEMPLATE_SWITCHING,
                ((first, InputRole.TEMPLATE), (tso, InputRole.STRUCTURAL),
                 (reverse_transcriptase, InputRole.CATALYST)), (switched,))

    preamp_polymerase = graph.input("preamp_polymerase", "full-length cDNA PCR polymerase",
                                    MaterialKind.CATALYST)
    ispcr = graph.input("ispcr", "ISPCR primer used at both cDNA ends",
                        MaterialKind.OLIGO)
    preamp_geometry = PCRGeometry.symmetric(
        preamp_template, PCRPrimer("ISPCR", rt.SMART_HANDLE),
        statement=("The source-defined SMART handle occurs in opposite orientations "
                   "at the two cDNA ends, so one ISPCR reagent supplies both inward sites."))
    preamp_inputs = ((switched, InputRole.TEMPLATE),
                     (ispcr, InputRole.STRUCTURAL),
                     (preamp_polymerase, InputRole.CATALYST))
    primer_reviewed = ["ispcr"]
    amplified = graph.product("amplified", "amplified full-length cDNA",
                              MaterialKind.MOLECULAR, state=preamp_template)
    graph.react("preamp", "PCR-amplify full-length cDNA", StepKind.PCR,
                preamp_inputs, (amplified,), geometry=preamp_geometry)

    transposome = graph.input("transposome", "s5/s7-loaded Nextera transposome",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagmented_product = graph.product(
        "tagmented", "heterologous s5/s7-tagmented cDNA fragments",
        MaterialKind.MOLECULAR, state=tagmented("s5", "s7"))
    tagment_waste = graph.product(
        "tagment_waste", "s5/s5, s7/s7 and other non-library Tn5 products",
        MaterialKind.WASTE)
    graph.react("tagment", "Tagment amplified cDNA with Nextera s5/s7 adapters",
                StepKind.TAGMENTATION,
                ((amplified, InputRole.TEMPLATE),
                 (transposome, InputRole.STRUCTURAL), (tn5, InputRole.CATALYST)),
                (tagmented_product, tagment_waste))
    graph.discard(tagment_waste,
                  "lacks the heterologous s5/s7 end combination required for amplification")

    gap_polymerase = graph.input("gap_polymerase", "Nextera nine-base gap-fill polymerase",
                                 MaterialKind.CATALYST)
    filled = graph.product("gap_filled", "gap-filled heterologous Nextera fragments",
                           MaterialKind.MOLECULAR, state=gap_filled())
    graph.react("gap_fill", "Fill both nine-base Tn5 gaps at 72 degrees C",
                StepKind.EXTENSION,
                ((tagmented_product, InputRole.TEMPLATE),
                 (gap_polymerase, InputRole.CATALYST)), (filled,))

    final = final_library()
    index_forward = graph.input("index_forward", "Nextera N5xx index primer",
                                MaterialKind.OLIGO)
    index_reverse = graph.input("index_reverse", "Nextera N7xx index primer",
                                MaterialKind.OLIGO)
    index_polymerase = graph.input("index_polymerase", "Nextera index-PCR polymerase",
                                   MaterialKind.CATALYST)
    index_geometry = PCRGeometry.exact(
        final, PCRPrimer("N5xx", nx.n5xx_primer("N" * 8)),
        PCRPrimer("N7xx", nx.n7xx_primer("N" * 8)),
        statement=("The canonical N5xx and N7xx primers face inward from the s5 and s7 "
                   "entry points and install the indexed P5/P7 ends."))
    primers = SEQ_PRIMERS
    main_library = graph.product("library", protocol + " sequencing-library molecules",
                                 MaterialKind.LIBRARY, state=final)
    graph.react("index_pcr", "Amplify and install indexed P5/P7 ends",
                StepKind.PCR,
                ((filled, InputRole.TEMPLATE),
                 (index_forward, InputRole.STRUCTURAL),
                 (index_reverse, InputRole.STRUCTURAL),
                 (index_polymerase, InputRole.CATALYST)),
                (main_library,), geometry=index_geometry)
    graph.sequence(main_library, SequencingPlan(
        library=final, primers=primers,
        acquired_roles=tuple(p.role for p in primers)))

    reviewed = [
        "mrna", "rt_primer", "reverse_transcriptase", "reverse_transcription",
        "tso", "template_switch", *primer_reviewed, "preamp_polymerase", "preamp",
        "transposome", "tn5", "tagment", "gap_polymerase", "gap_fill",
        "index_forward", "index_reverse", "index_polymerase", "index_pcr",
        "library",
    ]
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()

def preamp_template(insert_nt: int = 36, *, linker_inferred: bool = True) -> Construct:
    """Full-length cDNA before ISPCR.

    The ACAT linker is inferred when a descendant names Smart-seq2 without printing
    its TSO. Pass ``linker_inferred=False`` only when that source prints ACAT itself.
    """
    if insert_nt < 1:
        raise ValueError("SMART-seq2 cDNA insert must be non-empty")
    return Construct([
        Segment("5′ SMART handle", rt.SMART_HANDLE, "tso"),
        Segment("TSO linker", TSO_LINKER, inferred=linker_inferred),
        Segment("TSO G contribution", "GGG", "tso"),
        Segment("transcript cDNA", "X" * insert_nt, placeholder=True),
        Segment("poly(A) copy", "A" * DT_LEN),
        Segment("oligo-dT linker reverse complement", revcomp(DT_LINKER)),
        Segment("3′ SMART handle reverse complement", revcomp(rt.SMART_HANDLE), "tso"),
    ], name="SMART-seq2 full-length cDNA before ISPCR")


def ispcr_geometry(template: Construct | None = None) -> PCRGeometry:
    template = template or preamp_template()
    return PCRGeometry.symmetric(
        template, PCRPrimer("ISPCR", ISPCR),
        statement=("SMART-seq2 places the same 23-nt ISPCR handle at both cDNA ends, "
                   "so one primer species supplies both inward-facing sites."))


def nextera_index_geometry(final_library: Construct) -> PCRGeometry:
    """Invariant N5xx/N7xx landing sites; representative indices do not affect binding."""
    return PCRGeometry.exact(
        final_library,
        PCRPrimer("N5xx", nx.n5xx_primer("A" * 8), anneal_nt=len(nx.S5)),
        PCRPrimer("N7xx", nx.n7xx_primer("A" * 8), anneal_nt=len(nx.S7)),
        statement=("Nextera index primers terminate in s5 and s7; those 3′ regions face "
                   "inward on unlike-ended tagmented cDNA. Representative A8 index "
                   "windows show the invariant landing geometry."))
