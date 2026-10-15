"""Zilionis et al. inDrop v2 reaction path and source-owned oligo inputs."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parents[1] / "indrop-v1__10.1016+j.cell.2015.04.044" / "tools"),
                str(HERE.parents[1] / "lib")]

import illumina as il
import indrop as I
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

SOURCE = "Zilionis et al. 2017, doi:10.1038/nprot.2016.154; Supplementary Table 1"
PE2_N6 = I.PE2 + "N" * 6
PE2_PCR = il.P5 + "GGTC" + I.PE2[:18]
PE1_INDEX1 = il.P7 + "CGTGAT" + I.PE1[:16]


def protocol_graph() -> ProtocolGraph:
    """Trace v2 through bead synthesis, droplets, IVT, conversion and indexed PCR."""
    graph = ProtocolGraph(I.V2, authority=SOURCE)
    seed = graph.input("bead_seed", "acrydite photocleavable bead primer",
                       MaterialKind.OLIGO, state=I.SOURCE_ACRYDITE_BASES, source=SOURCE)
    plate = graph.input("barcode_plate_oligos", "two split-pool barcode plate oligo sets",
                        MaterialKind.OLIGO, source=SOURCE)
    bead_pol = graph.input("bead_polymerase", "bead-primer extension polymerase",
                           MaterialKind.CATALYST)
    bead = graph.product("barcoded_bead", "finished split-pool barcode bead primer",
                         MaterialKind.OLIGO, state=I.bead_primer(), source=SOURCE)
    graph.react("bead_synthesis", "Build the bead barcode in two split-pool extensions",
                StepKind.EXTENSION,
                ((seed, InputRole.TEMPLATE), (plate, InputRole.STRUCTURAL),
                 (bead_pol, InputRole.CATALYST)), (bead,))

    cells = graph.input("cells", "single-cell suspension", MaterialKind.BIOLOGICAL,
                        primary=True)
    loading = graph.input("droplet_loading", "one-cell/one-bead droplet loading condition",
                          MaterialKind.CONDITION)
    droplets = graph.product("productive_droplets", "one-cell/one-bead droplets",
                             MaterialKind.MOLECULAR)
    excluded = graph.product("excluded_droplets", "empty or multiplet droplets",
                             MaterialKind.WASTE)
    graph.react("partition", "Co-encapsulate cells and beads and photorelease primers",
                StepKind.PARTITION,
                ((cells, InputRole.TEMPLATE), (bead, InputRole.SPIKE_IN),
                 (loading, InputRole.CONDITION)), (droplets, excluded))
    graph.discard(excluded, "empty and multiplet droplets are not accepted cell libraries")
    rt = graph.input("reverse_transcriptase", "reverse transcriptase",
                     MaterialKind.CATALYST)
    first = graph.product("barcoded_cdna", "bead-barcoded first-strand cDNA",
                          MaterialKind.MOLECULAR)
    graph.react("reverse_transcription", "Reverse-transcribe from released T19V primers",
                StepKind.REVERSE_TRANSCRIPTION,
                ((droplets, InputRole.TEMPLATE), (bead, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (first,))
    second_pol = graph.input("second_strand_polymerase", "second-strand polymerase",
                             MaterialKind.CATALYST)
    ds = graph.product("t7_templates", "double-stranded T7-bearing cDNA",
                       MaterialKind.MOLECULAR)
    graph.react("second_strand", "Make the T7 promoter double-stranded",
                StepKind.EXTENSION,
                ((first, InputRole.TEMPLATE), (second_pol, InputRole.CATALYST)), (ds,))
    t7 = graph.input("t7_polymerase", "T7 RNA polymerase", MaterialKind.CATALYST)
    arna = graph.product("arna", "linearly amplified antisense RNA",
                         MaterialKind.MOLECULAR)
    graph.react("ivt", "Transcribe pooled cDNA from the T7 promoter",
                StepKind.TRANSCRIPTION,
                ((ds, InputRole.TEMPLATE), (t7, InputRole.CATALYST)), (arna,))

    pe2n6 = graph.input("pe2_n6", "PE2-N6 random primer", MaterialKind.OLIGO,
                        state=PE2_N6, source=SOURCE)
    converted = graph.product("sense_cdna", "PE2-tagged sense cDNA",
                              MaterialKind.MOLECULAR)
    graph.react("second_rt", "Random-prime aRNA with PE2-N6",
                StepKind.REVERSE_TRANSCRIPTION,
                ((arna, InputRole.TEMPLATE), (pe2n6, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (converted,))
    forward = graph.input("pcr_forward", "PE2-side library PCR primer",
                          MaterialKind.OLIGO, state=PE2_PCR, source=SOURCE)
    reverse = graph.input("pcr_reverse", "PE1 index-1 library PCR primer",
                          MaterialKind.OLIGO, state=PE1_INDEX1, source=SOURCE)
    pcr_pol = graph.input("pcr_polymerase", "library PCR polymerase",
                          MaterialKind.CATALYST)
    final = I.v2_library()
    library = graph.product("library", final.name, MaterialKind.LIBRARY, state=final)
    graph.react("library_pcr", "Add P5, P7 and sample index by enrichment PCR",
                StepKind.PCR,
                ((converted, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (pcr_pol, InputRole.CATALYST)),
                (library,), geometry=PCRGeometry.unavailable(
                    final, reason="The pre-PCR IVT-conversion intermediate is not represented base-for-base."))
    graph.sequence(library, SequencingPlan(
        library=final, primers=I.V2_PRIMERS,
        acquired_roles=tuple(p.role for p in I.V2_PRIMERS),
        read_lengths={"Index 1 (i7)": 6, "Read 2": 51},
        read_uses={"Read 1": "include", "Index 1 (i7)": "demultiplex",
                   "Read 2": "extract"},
        limitation="Read 1 length is run-dependent."))
    graph.declare_reviewed_scope((
        "bead_seed", "barcode_plate_oligos", "bead_polymerase", "bead_synthesis",
        "cells", "droplet_loading", "partition", "reverse_transcriptase",
        "reverse_transcription", "second_strand_polymerase", "second_strand",
        "t7_polymerase", "ivt", "pe2_n6", "second_rt", "pcr_forward",
        "pcr_reverse", "pcr_polymerase", "library_pcr", "library"), SOURCE)
    return graph.require_publishable()
