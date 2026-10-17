"""Sequence model for Smallwood et al. scBS-seq (Nature Methods, 2014).

The defining paper prints oligo1, oligo2 and PE1.0 and cites Quail et al. for the
indexed iPCRTag primer.  Quail's Supplementary Table 1 prints that primer family,
the PE adapter and the dedicated index-read primer.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Segment, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (Disposition, Terminal, InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)


def _seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


# Printed verbatim in the scBS-seq Online Methods.
OLIGO1_HANDLE = "CTACACGACGCTCTTCCGATCT"
OLIGO2_HANDLE = "TGCTGAACCGCTCTTCCGATCT"
RANDOM_NT = 9
PE1 = il.P5 + il.TRUSEQ_READ1[4:]  # the shared ACAC is present only once

# Quail et al. 2012, doi:10.1038/nmeth.1814, Supplementary Table 1.
IPCRTAG_INDEX_NT = 8
READ2_PRIMER = il.PE_READ2_PRIMER
IPCRTAG_3_ARM = "GAGAT" + READ2_PRIMER
IPCRTAG_T1 = il.P7 + "AACGTGAT" + IPCRTAG_3_ARM
INDEX_PRIMER = "AAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTC"


def oligo1_segments() -> list[Segment]:
    return [_seg("oligo1 handle", OLIGO1_HANDLE, "r1"),
            _seg("N9", "N" * RANDOM_NT, placeholder=True)]


def oligo2_segments() -> list[Segment]:
    return [_seg("oligo2 handle", OLIGO2_HANDLE, "r2"),
            _seg("N9", "N" * RANDOM_NT, placeholder=True)]


def ipcrtag_primer_segments() -> list[Segment]:
    """Quail iPCRTag primer family; the eight-base member index is variable."""
    return [
        _seg("P7", il.P7, "p7"),
        _seg("i7", "I" * IPCRTAG_INDEX_NT, "cbc", placeholder=True,
             feature=feature("cell_i7", "cell_barcode", "whitelist", whitelist="published iPCRTag primer set")),
        _seg("3' arm", IPCRTAG_3_ARM, "r2"),
    ]


def final_library(insert_nt: int = 28) -> Construct:
    """Indexed library, top strand 5'->3' from P5 to P7'.

    A short placeholder insert keeps the page readable.  The two N9s are random
    priming sequence, not UMIs, despite sharing the UMI colour token.
    """
    insert = "X" * insert_nt if insert_nt <= 20 else "XXXXXXXXXX...XXXXXXXXXX"
    lib = Construct([
        # P5 and Read 1 overlap by ACAC in PE1, so segment the 58-nt primer this way.
        _seg("P5", il.P5, "p5"),
        _seg("Read 1 remainder", il.TRUSEQ_READ1[4:], "r1"),
        _seg("oligo1 N9", "N" * RANDOM_NT, placeholder=True),
        _seg("bisulfite insert", insert, placeholder=True),
        _seg("oligo2 N9", "N" * RANDOM_NT, placeholder=True),
        _seg("index-read site", revcomp(IPCRTAG_3_ARM), "r2"),
        _seg("i7", "I" * IPCRTAG_INDEX_NT, "cbc", placeholder=True,
             feature=feature("cell_i7", "cell_barcode", "whitelist", whitelist="published iPCRTag primer set")),
        _seg("P7'", il.P7_RC, "p7"),
    ], name="scBS-seq indexed library")
    problems = sp.verify(lib, SEQ_PRIMERS,
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if problems:
        raise ValueError("invalid scBS-seq final library: " + "; ".join(problems))
    return lib


SEQ_PRIMERS = [
    sp.custom("Read 1", "TruSeq Read 1", il.TRUSEQ_READ1,
              'Illumina "Illumina Adapter Sequences" #1000000002694'),
    sp.custom("Index 1 (i7)", "historical iPCRTag index primer", INDEX_PRIMER,
              "Quail et al. 2012 Supplementary Table 1"),
    sp.custom("Read 2", "historical paired-end Read 2", READ2_PRIMER,
              'Illumina "Oligonucleotide Sequences for Paired End DNA" (obsolete)'),
]


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1038/nmeth.3035"
    graph = ProtocolGraph("scBS-seq", authority=source)
    genome = graph.input("genome", "single-cell genomic DNA",
                         MaterialKind.BIOLOGICAL, primary=True)
    bisulfite = graph.input("bisulfite", "bisulfite-conversion conditions",
                            MaterialKind.CONDITION)
    converted = graph.product("converted_dna", "fragmented bisulfite-converted DNA",
                              MaterialKind.MOLECULAR)
    graph.react("bisulfite_conversion", "Bisulfite-convert and fragment genomic DNA",
                StepKind.CHEMICAL_CONVERSION,
                ((genome, InputRole.TEMPLATE), (bisulfite, InputRole.CONDITION)),
                (converted,))
    oligo1 = graph.input("oligo1", "biotinylated oligo1 random primer",
                         MaterialKind.OLIGO,
                         state=OLIGO1_HANDLE + "N" * RANDOM_NT,
                         source=source + "; Online Methods, [Btn]oligo1")
    klenow = graph.input("klenow", "Klenow exo- polymerase",
                         MaterialKind.CATALYST)
    first = graph.product("oligo1_products", "displaced oligo1 extension products",
                          MaterialKind.MOLECULAR)
    graph.react("oligo1_extension", "Random-prime for five extension cycles",
                StepKind.EXTENSION,
                ((converted, InputRole.TEMPLATE), (oligo1, InputRole.STRUCTURAL),
                 (klenow, InputRole.CATALYST)), (first,))
    exo = graph.input("exonuclease_i", "Exonuclease I for unused oligo1",
                      MaterialKind.CATALYST)
    cleaned = graph.product("primer_cleaned", "oligo1 products after free-primer digestion",
                            MaterialKind.MOLECULAR)
    primer_waste = graph.product("free_oligo1_waste", "digested free oligo1",
                                 MaterialKind.MOLECULAR)
    graph.react("remove_free_oligo1", "Digest unused oligo1 with Exonuclease I",
                StepKind.DIGESTION,
                ((first, InputRole.TEMPLATE), (exo, InputRole.CATALYST)),
                (cleaned, primer_waste))
    graph.mark_terminal(Terminal(primer_waste.id, Disposition.UNRESOLVED, "Digested primer products are not separately tracked through the subsequent AMPure cleanup."))
    ampure = graph.input("ampure", "AMPure XP cleanup beads", MaterialKind.SELECTION)
    purified = graph.product("purified_oligo1_products", "purified oligo1 extension products",
                             MaterialKind.MOLECULAR)
    ampure_waste = graph.product("ampure_waste", "removed short fragments and buffer",
                                 MaterialKind.WASTE)
    graph.react("ampure_cleanup", "Purify oligo1 products with AMPure XP",
                StepKind.CLEANUP,
                ((cleaned, InputRole.TEMPLATE), (ampure, InputRole.SELECTION)),
                (purified, ampure_waste))
    graph.discard(ampure_waste, "unretained material is removed during AMPure purification")
    beads = graph.input("streptavidin_beads", "streptavidin capture beads",
                        MaterialKind.SELECTION)
    captured = graph.product("captured_products", "bead-bound oligo1 products",
                             MaterialKind.MOLECULAR)
    rejected = graph.product("uncaptured_material", "uncaptured conversion products",
                             MaterialKind.WASTE)
    graph.react("capture", "Capture biotinylated products and wash",
                StepKind.AFFINITY_CAPTURE,
                ((purified, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (captured, rejected))
    graph.discard(rejected, "unbound material is removed during streptavidin washes")
    oligo2 = graph.input("oligo2", "oligo2 random primer", MaterialKind.OLIGO,
                         state=OLIGO2_HANDLE + "N" * RANDOM_NT,
                         source=source + "; Online Methods, oligo2")
    second = graph.product("dual_handle_products", "dual-handle extension products",
                           MaterialKind.MOLECULAR)
    graph.react("oligo2_extension", "Random-prime the captured strand with oligo2",
                StepKind.EXTENSION,
                ((captured, InputRole.TEMPLATE), (oligo2, InputRole.STRUCTURAL),
                 (klenow, InputRole.CATALYST)), (second,))
    forward = graph.input("pe1_primer", "PE1.0 forward PCR primer", MaterialKind.OLIGO,
                          state=PE1, source=source + "; Online Methods, PE1.0")
    reverse = graph.input("ipcrtag_primer", "indexed iPCRTag reverse primer",
                          MaterialKind.OLIGO, state=IPCRTAG_T1,
                          source="Quail et al. 2012, doi:10.1038/nmeth.1814, Supplementary Table 1 iPCRtagT1")
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    library = graph.product("library", final_library().name, MaterialKind.LIBRARY,
                            state=final_library())
    graph.react("indexed_pcr", "Add flow-cell arms and the cell index: 12–13 PCR cycles",
                StepKind.PCR,
                ((second, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final_library(), reason="The intermediate random-primer product is represented schematically; the published primers and final arms are exact."))
    graph.sequence(library, SequencingPlan(
        library=final_library(), primers=tuple(SEQ_PRIMERS),
        acquired_roles=("Read 1", "Index 1 (i7)", "Read 2"),
        read_lengths={"Read 1": 100, "Index 1 (i7)": IPCRTAG_INDEX_NT,
                      "Read 2": 100},
        read_uses={"Read 1": "include", "Index 1 (i7)": "demultiplex",
                   "Read 2": "include"}))
    graph.declare_reviewed_scope((
        "genome", "bisulfite", "bisulfite_conversion", "oligo1", "klenow",
        "oligo1_extension", "exonuclease_i", "remove_free_oligo1",
        "ampure", "ampure_cleanup", "streptavidin_beads", "capture", "oligo2",
        "oligo2_extension", "pe1_primer", "ipcrtag_primer", "pcr_polymerase",
        "indexed_pcr", "library"), source)
    return graph.require_publishable()
