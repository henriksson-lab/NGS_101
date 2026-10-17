"""Constructors for scNOMe-seq and scCOOL-seq library molecules."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Segment, feature, revcomp
import illumina as il
import seqprimers as sp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

SCNOME = "scNOMe-seq"
SCCOOL = "scCOOL-seq"

RANDOM_NT = 9
COOL_P1_HANDLE = il.TRUSEQ_READ1[-22:]
COOL_P2_HANDLE = il.TRUSEQ_READ2[-21:]
FORWARD_P5 = il.P5 + il.TRUSEQ_READ1[4:]
INDEX_NT_NOME = 6
INDEX_NT_COOL = 8       # standard indexed-primer-shaped model; exact NEB set not printed

SOURCE_FORWARD_P5 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
SOURCE_COOL_P1 = "CTACACGACGCTCTTCCGATCTNNNNNNNNN"


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def primer1_segments() -> list[Segment]:
    return [seg("Read 1 handle", COOL_P1_HANDLE, "r1"),
            seg("N9", "N" * RANDOM_NT, placeholder=True)]


def primer2_segments() -> list[Segment]:
    return [seg("Read 2 handle", COOL_P2_HANDLE, "r2"),
            seg("N9", "N" * RANDOM_NT, placeholder=True)]


SEQ_PRIMERS = (
    sp.custom("Read 1", "TruSeq Read 1", il.TRUSEQ_READ1,
              "Illumina adapter sequence"),
    sp.custom("Index 1 (i7)", "TruSeq Index 1", il.INDEX1_PRIMER,
              "Illumina adapter sequence"),
    sp.custom("Read 2", "TruSeq Read 2", il.TRUSEQ_READ2,
              "Illumina adapter sequence"),
)
SCNOME_READ_LENGTHS = {"Read 1": 100, "Index 1 (i7)": 6, "Read 2": 100}
SCNOME_READ_USES = {"Read 1": "include", "Read 2": "include"}


def nome_library() -> Construct:
    """Outer published primer arms with one explicit unknown kit-derived junction each."""
    con = Construct([
        seg("P5", il.P5, "p5"), seg("Read 1 remainder", il.TRUSEQ_READ1[4:], "r1"),
        seg("kit remnant", "?", inferred=True, placeholder=True,
            note="Pico Methyl-Seq internal primer sequence and length unresolved in this source pass"),
        seg("bisulfite insert", "XXXXXXXX...XXXXXXXX", placeholder=True),
        seg("kit remnant'", "?", inferred=True, placeholder=True,
            note="Pico Methyl-Seq internal primer sequence and length unresolved in this source pass"),
        seg("Read 2 arm'", revcomp(il.TRUSEQ_READ2), "r2"),
        seg("i7'", "I" * INDEX_NT_NOME, "cbc", placeholder=True,
            feature=feature("sample_i7", "sample_index", "unknown")),
        seg("P7'", il.P7_RC, "p7"),
    ], name="scNOMe-seq library")
    _validate(con)
    return con


def cool_library() -> Construct:
    arm = revcomp(il.TRUSEQ_READ2)
    known = revcomp(COOL_P2_HANDLE)
    if not arm.startswith(known):
        raise ValueError("primer-2 handle is not a suffix of TruSeq Read 2")
    con = Construct([
        seg("P5", il.P5, "p5", inferred=True),
        seg("Read 1 remainder", il.TRUSEQ_READ1[4:], "r1", inferred=True),
        seg("primer-1 N9", "N" * RANDOM_NT, placeholder=True),
        seg("bisulfite insert", "XXXXXXXX...XXXXXXXX", placeholder=True),
        seg("primer-2 N9'", "N" * RANDOM_NT, placeholder=True),
        seg("primer-2 handle'", known, "r2"),
        seg("reverse-primer extension", arm[len(known):], "r2", inferred=True,
            note="standard TruSeq-shaped NEB indexed primer; exact product not printed"),
        seg("i7'", "I" * INDEX_NT_COOL, "cbc", placeholder=True, inferred=True,
            feature=feature("sample_i7", "sample_index", "unknown")),
        seg("P7'", il.P7_RC, "p7", inferred=True),
    ], name="scCOOL-seq library")
    _validate(con)
    return con


def _validate(con: Construct) -> None:
    for primer in SEQ_PRIMERS:
        if sp.locate(con, primer) is None:
            raise ValueError(f"{con.name}: {primer.name} site is absent")


def library(protocol: str) -> Construct:
    if protocol == SCNOME: return nome_library()
    if protocol == SCCOOL: return cool_library()
    raise ValueError(protocol)


def protocol_graph(protocol: str) -> ProtocolGraph:
    doi = "10.7554/eLife.23203" if protocol == SCNOME else "10.1038/cr.2017.82"
    source = f"Defining source: https://doi.org/{doi}"
    graph = ProtocolGraph(protocol, authority=source)
    genome = graph.input("genome", "single-cell genomic DNA",
                         MaterialKind.BIOLOGICAL, primary=True)
    mcvi = graph.input("mcvi", "M.CviPI GpC methyltransferase",
                       MaterialKind.CATALYST)
    marked = graph.product("gpc_marked_dna", "accessibility-marked genomic DNA",
                           MaterialKind.MOLECULAR)
    graph.react("gpc_marking", "Methylate accessible GpC sites",
                StepKind.CHEMICAL_CONVERSION,
                ((genome, InputRole.TEMPLATE), (mcvi, InputRole.CATALYST)), (marked,))
    bisulfite = graph.input("bisulfite", "bisulfite-conversion conditions",
                            MaterialKind.CONDITION)
    converted = graph.product("converted_dna", "bisulfite-converted marked DNA",
                              MaterialKind.MOLECULAR)
    graph.react("bisulfite_conversion", "Bisulfite-convert marked DNA",
                StepKind.CHEMICAL_CONVERSION,
                ((marked, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)),
                (converted,))

    if protocol == SCNOME:
        kit = graph.input("pico_kit_primers", "Pico Methyl-Seq internal primers",
                          MaterialKind.OLIGO, inferred=True,
                          source=source + " (proprietary internal oligos unavailable)")
        polymerase = graph.input("kit_polymerase", "Pico Methyl-Seq polymerase",
                                 MaterialKind.CATALYST)
        preamp = graph.product("pico_products", "Pico Methyl-Seq preamplified products",
                               MaterialKind.MOLECULAR, inferred=True,
                               source=source + " (internal junction bases unavailable)")
        graph.react("pico_library_generation", "Pico Methyl-Seq low-input preparation: 1:2 random-primer dilution, 10 first-amplification cycles",
                    StepKind.EXTENSION,
                    ((converted, InputRole.TEMPLATE), (kit, InputRole.STRUCTURAL),
                     (polymerase, InputRole.CATALYST)), (preamp,))
        forward = graph.input("forward_p5", "published P5 forward primer",
                              MaterialKind.OLIGO,
                              state=Construct([seg("forward_p5", FORWARD_P5)]))
        reverse = graph.input("indexed_p7", "published indexed P7 primer",
                              MaterialKind.OLIGO)
        pcr_template = preamp
        pcr_id = "indexed_pcr"
        pcr_label = "Add published P5 and indexed P7 arms by PCR"
        reviewed = ("genome", "mcvi", "gpc_marking", "bisulfite",
                    "bisulfite_conversion", "pico_kit_primers", "kit_polymerase",
                    "pico_library_generation", "forward_p5", "indexed_p7")
    else:
        p1 = graph.input("random_primer1", "biotinylated random primer 1",
                         MaterialKind.OLIGO, state=tuple(primer1_segments()))
        polymerase = graph.input("extension_polymerase", "random-priming polymerase",
                                 MaterialKind.CATALYST)
        first = graph.product("first_strands", "biotinylated first-strand products",
                              MaterialKind.MOLECULAR)
        graph.react("first_priming", "Random-prime the first strand",
                    StepKind.EXTENSION,
                    ((converted, InputRole.TEMPLATE), (p1, InputRole.STRUCTURAL),
                     (polymerase, InputRole.CATALYST)), (first,))
        beads = graph.input("streptavidin_beads", "streptavidin capture beads",
                            MaterialKind.SELECTION)
        captured = graph.product("captured_strands", "bead-bound first strands",
                                 MaterialKind.MOLECULAR)
        waste = graph.product("uncaptured_material", "uncaptured material",
                              MaterialKind.WASTE)
        graph.react("capture", "Capture biotinylated first strands and wash",
                    StepKind.AFFINITY_CAPTURE,
                    ((first, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                    (captured, waste))
        graph.discard(waste, "unbound material and removed converted templates are washed away")
        p2 = graph.input("random_primer2", "random primer 2", MaterialKind.OLIGO,
                         state=tuple(primer2_segments()))
        preamp = graph.product("dual_handle_products", "dual-handle second-strand products",
                               MaterialKind.MOLECULAR)
        graph.react("second_priming", "Random-prime the captured second strand",
                    StepKind.EXTENSION,
                    ((captured, InputRole.TEMPLATE), (p2, InputRole.STRUCTURAL),
                     (polymerase, InputRole.CATALYST)), (preamp,))
        forward = graph.input("forward_pe1", "Forward PE1.0 PCR primer",
                              MaterialKind.OLIGO, inferred=True,
                              source=source + " (sequence inferred from reagent name)")
        reverse = graph.input("indexed_reverse", "NEB indexed reverse PCR primer",
                              MaterialKind.OLIGO, inferred=True,
                              source=source + " (exact primer sequence unavailable)")
        pcr_template = preamp
        pcr_id = "indexed_pcr"
        pcr_label = "Add flow-cell arms and index by PCR"
        reviewed = ("genome", "mcvi", "gpc_marking", "bisulfite",
                    "bisulfite_conversion", "random_primer1", "extension_polymerase",
                    "first_priming", "streptavidin_beads", "capture", "random_primer2",
                    "second_priming", "forward_pe1", "indexed_reverse")
    pcr_pol = graph.input("pcr_polymerase", "library PCR polymerase",
                          MaterialKind.CATALYST)
    final = library(protocol)
    # The scNOMe construct deliberately contains literal question-mark junctions; keep
    # those in the dedicated final panel rather than asking the graph renderer to invent
    # a complementary base for them.
    graph_state = None if protocol == SCNOME else final
    out = graph.product("library", final.name,
                        MaterialKind.MOLECULAR if protocol == SCNOME else MaterialKind.LIBRARY, state=graph_state,
                        inferred=protocol == SCCOOL,
                        source=(source + " (outer PCR arms inferred from named primers)"
                                if protocol == SCCOOL else source))
    graph.react(pcr_id, pcr_label, StepKind.PCR,
                ((pcr_template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (pcr_pol, InputRole.CATALYST)),
                (out,), geometry=PCRGeometry.unavailable(
                    final, reason="The low-input intermediate does not disclose enough sequence to locate both PCR sites exactly."))
    extra_reviewed = ()
    if protocol == SCNOME:
        ampure = graph.input("ampure", "Agencourt AMPureXP beads, 1:1 beads:library ratio",
                             MaterialKind.SELECTION)
        cleaned = graph.product("cleaned_library", "AMPure-purified scNOMe-seq library",
                                MaterialKind.LIBRARY)
        waste = graph.product("ampure_waste", "unrecovered material from AMPure library purification",
                              MaterialKind.WASTE)
        graph.react("ampure_cleanup", "Purify amplified libraries with AMPureXP before Bioanalyzer QC and pooling",
                    StepKind.CLEANUP,
                    ((out, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)), (cleaned, waste))
        graph.discard(waste, "not recovered in the library eluate; no exact composition assigned")
        out = cleaned
        extra_reviewed = ("ampure", "ampure_cleanup", "cleaned_library")
    graph.sequence(out, SequencingPlan(
        library=final, primers=tuple(SEQ_PRIMERS),
        acquired_roles=("Read 1", "Index 1 (i7)", "Read 2"),
        read_lengths=(SCNOME_READ_LENGTHS if protocol == SCNOME else
                      {"Read 1": 150, "Read 2": 150}),
        read_uses=(SCNOME_READ_USES if protocol == SCNOME else
                   {"Read 1": "include", "Index 1 (i7)": "demultiplex",
                    "Read 2": "include"}),
        limitation="Internal commercial-kit junction sequences are unavailable."))
    graph.declare_reviewed_scope((*reviewed, *extra_reviewed, "pcr_polymerase", pcr_id, "library"), source)
    return graph.require_publishable()


def mark_and_convert_graph(graph: ProtocolGraph):
    """The shared GpC marking and bisulfite conversion front end."""
    genome = graph.input("genome", "single-cell genomic DNA", MaterialKind.BIOLOGICAL,
                         primary=True)
    mcvi = graph.input("mcvi", "M.CviPI GpC methyltransferase", MaterialKind.CATALYST)
    marked = graph.product("gpc_marked_dna", "accessibility-marked genomic DNA",
                           MaterialKind.MOLECULAR)
    graph.react("gpc_marking", "Methylate accessible GpC sites", StepKind.CHEMICAL_CONVERSION,
                ((genome, InputRole.TEMPLATE), (mcvi, InputRole.CATALYST)), (marked,))
    bisulfite = graph.input("bisulfite", "bisulfite-conversion conditions",
                            MaterialKind.CONDITION)
    converted = graph.product("converted_dna", "bisulfite-converted marked DNA",
                              MaterialKind.MOLECULAR)
    graph.react("bisulfite_conversion", "Bisulfite-convert marked DNA",
                StepKind.CHEMICAL_CONVERSION,
                ((marked, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)),
                (converted,))
    return converted


def finish_indexed_graph(graph: ProtocolGraph, template, forward, reverse,
                         *, source: str, reviewed: tuple[str, ...]):
    """Shared indexed PCR/readout boundary after a protocol's priming sequence."""
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    final = library(SCCOOL)
    out = graph.product("library", final.name, MaterialKind.LIBRARY, state=final,
                        inferred=True,
                        source=source + " (outer PCR arms inferred from named primers)")
    graph.react("indexed_pcr", "Add flow-cell arms and index by PCR", StepKind.PCR,
                ((template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (out,), geometry=PCRGeometry.unavailable(
                    final, reason="The low-input intermediate does not disclose enough sequence to locate both PCR sites exactly."))
    graph.sequence(out, SequencingPlan(
        library=final, primers=tuple(SEQ_PRIMERS),
        acquired_roles=("Read 1", "Index 1 (i7)", "Read 2"),
        read_lengths={"Read 1": 150, "Read 2": 150},
        read_uses={"Read 1": "include", "Index 1 (i7)": "demultiplex", "Read 2": "include"},
        limitation="Internal commercial-kit junction sequences are unavailable."))
    graph.declare_reviewed_scope((*reviewed, "pcr_polymerase", "indexed_pcr", "library"), source)
    return graph.require_publishable()
