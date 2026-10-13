"""Molecular constructors for the published inDrop v1 and v2 protocols."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Segment, feature, revcomp
import illumina as il
import seqprimers as sp
from adapter_ligation import AdapterLigationGeometry
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

V1 = "inDrop v1"
V2 = "inDrop v2"

LEADER = "CGATGACG"
T7_PROMOTER = "TAATACGACTCACTATAGGG"
MID = "ATACCACCATGG"
PE1 = il.TRUSEQ_READ1[3:]
W1_STAR = "AAGGCGTCACAAGCAATCACTC"
W1 = revcomp(W1_STAR)
BC2_LEN = 8
UMI_LEN = 6
DT_LEN = 19
PE2 = "TCGGCATTCCTGCTGAACCGCTCTTCCGATCT"

RLO = revcomp(PE2[2:])
V1_R2 = "CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"
V1_R1 = il.TRUSEQ_READ1[4:]
V2_R1 = PE2[2:]
V2_INDEX = revcomp(PE1)
V2_R2 = PE1

SOURCE_ACRYDITE_BASES = "CGATGACGTAATACGACTCACTATAGGGATACCACCATGGCTCTTTCCCTACACGACGCTCTTC"
SOURCE_PE2_N6 = "TCGGCATTCCTGCTGAACCGCTCTTCCGATCTNNNNNN"
BC1_FEATURE = feature("cell_barcode_1", "cell_barcode", "combinatorial",
                      group="cell_barcode", part="barcode 1")
BC2_FEATURE = feature("cell_barcode_2", "cell_barcode", "combinatorial",
                      group="cell_barcode", part="barcode 2")
UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "unknown")
FEATURES = {"barcode 1": BC1_FEATURE, "barcode 1'": BC1_FEATURE,
            "barcode 2": BC2_FEATURE, "barcode 2'": BC2_FEATURE,
            "UMI": UMI_FEATURE, "i7 read": I7_FEATURE}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


def var(name: str, n: int, letter: str, tag: str | None = None,
        inferred: bool = False) -> Segment:
    return seg(name, letter * n, tag, placeholder=True, inferred=inferred)


def bead_primer(bc1_len: int = 8, *, inferred: bool = False) -> Construct:
    if bc1_len not in (8, 9, 10, 11):
        raise ValueError("inDrop barcode 1 must be 8-11 nt")
    return Construct([
        seg("leader", LEADER, inferred=inferred), seg("T7 promoter", T7_PROMOTER, "tso", inferred=inferred),
        seg("middle", MID, inferred=inferred), seg("PE1", PE1, "r2", inferred=inferred),
        var("barcode 1'", bc1_len, "A", "cbc", inferred), seg("W1", W1, "r3", inferred=inferred),
        var("barcode 2'", BC2_LEN, "B", "cbc", inferred), var("UMI", UMI_LEN, "N", "umi", inferred),
        seg("T19", "T" * DT_LEN, inferred=inferred), var("V", 1, "V", None, inferred),
    ], name="finished inDrop bead primer")


V1_PRIMERS = (
    sp.custom("Read 1", "v1 Read 1", V1_R1, "upstream secondary source"),
    sp.custom("Read 2", "v1 Read 2", V1_R2, "upstream secondary source"),
)
V2_PRIMERS = (
    sp.custom("Read 1", "Custom Read 1", V2_R1, "Nature Protocols supplementary oligo table"),
    sp.custom("Index 1 (i7)", "Custom Index Read", V2_INDEX, "Nature Protocols supplementary oligo table"),
    sp.custom("Read 2", "Custom Read 2", V2_R2, "Nature Protocols supplementary oligo table"),
)
V1_READ_LENGTHS = {"Read 1": 51}
V1_READ_USES = {"Read 1": "extract", "Read 2": "include"}


def v1_library(bc1_len: int = 8) -> Construct:
    """Entire v1 sequence model is secondary-source and therefore inseparably inferred."""
    inf = True
    con = Construct([
        seg("P5", il.P5, "p5", inferred=inf),
        seg("Read 1 site", V1_R1, "r1", inferred=inf),
        var("barcode 1'", bc1_len, "A", "cbc", inf), seg("W1", W1, "r3", inferred=inf),
        var("barcode 2'", BC2_LEN, "B", "cbc", inf), var("UMI", UMI_LEN, "N", "umi", inf),
        seg("T19", "T" * DT_LEN, inferred=inf), var("V", 1, "V", None, inf),
        var("antisense mRNA", 19, "X", "r1", inf), seg("RLO", RLO, "r2", inferred=inf),
        seg("RT2 extension", "GAGACCG", "r2", inferred=inf), seg("P7'", il.P7_RC, "p7", inferred=inf),
    ], name="inDrop v1 library")
    _validate(con, V1_PRIMERS)
    return con


def v2_library(bc1_len: int = 8) -> Construct:
    if bc1_len not in (8, 9, 10, 11): raise ValueError("barcode 1 must be 8-11 nt")
    con = Construct([
        seg("P5", il.P5, "p5"), seg("GGTC", "GGTC", "r1"), seg("PE2", PE2, "r1"),
        var("sense mRNA", 19, "X", "r1"), seg("A19", "A" * DT_LEN), var("B", 1, "B"),
        var("UMI", UMI_LEN, "N", "umi"), var("barcode 2", BC2_LEN, "B", "cbc"),
        seg("W1*", W1_STAR, "r3"), var("barcode 1", bc1_len, "A", "cbc"),
        seg("PE1*", revcomp(PE1), "r2"), var("i7 read", 6, "I", "cbc"),
        seg("P7'", il.P7_RC, "p7"),
    ], name="inDrop v2 library")
    _validate(con, V2_PRIMERS)
    return con


def _validate(con: Construct, primers) -> None:
    for p in primers:
        if sp.locate(con, p) is None:
            raise ValueError(f"{con.name}: {p.name} site missing")


def library(protocol: str) -> Construct:
    if protocol == V1: return v1_library()
    if protocol == V2: return v2_library()
    raise ValueError(protocol)


def protocol_graph(protocol: str) -> ProtocolGraph:
    doi = "10.1016/j.cell.2015.04.044" if protocol == V1 else "10.1038/nprot.2016.154"
    source = f"Defining source: https://doi.org/{doi}"
    if protocol == V1:
        boundary = source + " (v1 oligo bases follow the cited secondary reconstruction)"
    else:
        boundary = source
    graph = ProtocolGraph(protocol, authority=source)
    bead_seed = graph.input("bead_seed", "acrydite photocleavable bead primer",
                            MaterialKind.OLIGO, inferred=protocol == V1, source=boundary)
    plate_oligos = graph.input("barcode_plate_oligos", "two split-pool barcode plate oligo sets",
                               MaterialKind.OLIGO, inferred=protocol == V1, source=boundary)
    bead_polymerase = graph.input("bead_polymerase", "bead-primer extension polymerase",
                                  MaterialKind.CATALYST)
    bead = graph.product("barcoded_bead", "finished split-pool barcode bead primer",
                         MaterialKind.OLIGO, state=bead_primer(inferred=protocol == V1),
                         inferred=protocol == V1, source=boundary)
    graph.react("bead_synthesis", "Build the bead barcode in two split-pool extensions",
                StepKind.EXTENSION,
                ((bead_seed, InputRole.TEMPLATE), (plate_oligos, InputRole.STRUCTURAL),
                 (bead_polymerase, InputRole.CATALYST)), (bead,))
    cells = graph.input("cells", "single-cell suspension", MaterialKind.BIOLOGICAL,
                        primary=True)
    loading = graph.input("droplet_loading", "one-cell/one-bead droplet loading condition",
                          MaterialKind.CONDITION)
    droplets = graph.product("productive_droplets", "droplets containing one cell and one barcode bead",
                             MaterialKind.MOLECULAR)
    excluded = graph.product("excluded_droplets", "empty or multiplet droplets", MaterialKind.WASTE)
    graph.react("partition", "Co-encapsulate cells and beads and photorelease primers",
                StepKind.PARTITION,
                ((cells, InputRole.TEMPLATE), (bead, InputRole.SPIKE_IN),
                 (loading, InputRole.CONDITION)), (droplets, excluded))
    graph.discard(excluded, "empty and multiplet droplets do not yield accepted single-cell libraries")
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    first = graph.product("barcoded_cdna", "bead-barcoded first-strand cDNA",
                          MaterialKind.MOLECULAR)
    graph.react("reverse_transcription", "Reverse-transcribe mRNA from released T19V primers",
                StepKind.REVERSE_TRANSCRIPTION,
                ((droplets, InputRole.TEMPLATE), (bead, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (first,))
    second_pol = graph.input("second_strand_polymerase", "second-strand polymerase",
                             MaterialKind.CATALYST)
    ds = graph.product("t7_templates", "double-stranded T7-bearing cDNA templates",
                       MaterialKind.MOLECULAR)
    graph.react("second_strand", "Make the bead-borne T7 promoter double-stranded",
                StepKind.EXTENSION,
                ((first, InputRole.TEMPLATE), (second_pol, InputRole.CATALYST)), (ds,))
    t7 = graph.input("t7_polymerase", "T7 RNA polymerase", MaterialKind.CATALYST)
    arna = graph.product("arna", "linearly amplified antisense RNA", MaterialKind.MOLECULAR)
    graph.react("ivt", "Linearly amplify pooled cDNA by T7 transcription",
                StepKind.TRANSCRIPTION,
                ((ds, InputRole.TEMPLATE), (t7, InputRole.CATALYST)), (arna,))

    if protocol == V1:
        fragmentation = graph.input("fragmentation_condition", "aRNA fragmentation condition",
                                    MaterialKind.CONDITION)
        fragmented = graph.product("fragmented_arna", "fragmented antisense RNA",
                                   MaterialKind.MOLECULAR)
        graph.react("fragment", "Fragment the amplified antisense RNA",
                    StepKind.FRAGMENTATION,
                    ((arna, InputRole.TEMPLATE), (fragmentation, InputRole.CONDITION)),
                    (fragmented,))
        rlo = graph.input("rlo", "blocked RNA ligation oligo", MaterialKind.ADAPTER,
                          inferred=True, source=boundary)
        ligase = graph.input("rna_ligase", "RNA ligase", MaterialKind.CATALYST)
        ligated = graph.product("rlo_ligated_arna", "RLO-ligated aRNA fragments",
                                MaterialKind.MOLECULAR, inferred=True, source=boundary)
        graph.react("rlo_ligation", "Ligate the blocked RLO to fragmented aRNA",
                    StepKind.LIGATION,
                    ((fragmented, InputRole.TEMPLATE), (rlo, InputRole.STRUCTURAL),
                     (ligase, InputRole.CATALYST)), (ligated,),
                    geometry=AdapterLigationGeometry.unavailable(
                        "The v1 RLO sequence and junction are available only from the secondary reconstruction."))
        conversion_primer = graph.input("second_rt_primer", "second reverse-transcription primer",
                                        MaterialKind.OLIGO, inferred=True, source=boundary)
        converted = graph.product("second_cdna", "second-round cDNA", MaterialKind.MOLECULAR,
                                  inferred=True, source=boundary)
        graph.react("second_rt", "Reverse-transcribe from the ligated RLO",
                    StepKind.REVERSE_TRANSCRIPTION,
                    ((ligated, InputRole.TEMPLATE), (conversion_primer, InputRole.STRUCTURAL),
                     (rt, InputRole.CATALYST)), (converted,))
        pcr_template = converted
        reviewed_tail = ("fragmentation_condition", "fragment", "rlo", "rna_ligase",
                         "rlo_ligation", "second_rt_primer", "second_rt")
        primers = V1_PRIMERS
        lengths = V1_READ_LENGTHS
        uses = V1_READ_USES
    else:
        pe2n6 = graph.input("pe2_n6", "PE2-N6 random primer", MaterialKind.OLIGO)
        converted = graph.product("sense_cdna", "PE2-tagged sense cDNA", MaterialKind.MOLECULAR)
        graph.react("second_rt", "Random-prime aRNA with PE2-N6",
                    StepKind.REVERSE_TRANSCRIPTION,
                    ((arna, InputRole.TEMPLATE), (pe2n6, InputRole.STRUCTURAL),
                     (rt, InputRole.CATALYST)), (converted,))
        pcr_template = converted
        reviewed_tail = ("pe2_n6", "second_rt")
        primers = V2_PRIMERS
        lengths = {"Index 1 (i7)": 6, "Read 2": 51}
        uses = {"Read 1": "include", "Index 1 (i7)": "demultiplex",
                "Read 2": "extract"}
    forward = graph.input("pcr_forward", "P5 library PCR primer", MaterialKind.OLIGO,
                          inferred=protocol == V1, source=boundary)
    reverse = graph.input("pcr_reverse", "indexed P7 library PCR primer", MaterialKind.OLIGO,
                          inferred=protocol == V1, source=boundary)
    pcr_pol = graph.input("pcr_polymerase", "library PCR polymerase", MaterialKind.CATALYST)
    final = library(protocol)
    out = graph.product("library", final.name, MaterialKind.LIBRARY, state=final,
                        inferred=protocol == V1, source=boundary)
    graph.react("library_pcr", "Add P5, P7 and sample index by enrichment PCR",
                StepKind.PCR,
                ((pcr_template, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (pcr_pol, InputRole.CATALYST)),
                (out,), geometry=PCRGeometry.unavailable(
                    final, reason="The pre-PCR IVT-conversion intermediate is not represented base-for-base."))
    graph.sequence(out, SequencingPlan(
        library=final, primers=tuple(primers),
        acquired_roles=tuple(p.role for p in primers), read_lengths=lengths,
        read_uses=uses,
        limitation=("All v1 final-library bases are secondary-source reconstructions."
                    if protocol == V1 else "Read 1 length is run-dependent.")))
    graph.declare_reviewed_scope((
        "bead_seed", "barcode_plate_oligos", "bead_polymerase", "bead_synthesis",
        "cells", "droplet_loading", "partition", "reverse_transcriptase",
        "reverse_transcription", "second_strand_polymerase", "second_strand",
        "t7_polymerase", "ivt", *reviewed_tail, "pcr_forward", "pcr_reverse",
        "pcr_polymerase", "library_pcr", "library"), source)
    return graph.require_publishable()

