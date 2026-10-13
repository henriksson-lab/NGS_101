"""Molecular construct model for FIPRESCI (Luo et al., 2023)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

ROUND1_FEATURE = feature("cell_barcode_round1", "cell_barcode", "combinatorial",
                         group="cell_barcode", part="indexed Tn5 round")
DROPLET_FEATURE = feature("cell_barcode_droplet", "cell_barcode", "combinatorial",
                          group="cell_barcode", part="droplet")
UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "fixed")
FEATURES = {"round-1 barcode": ROUND1_FEATURE,
            "round-1 barcode reverse complement": ROUND1_FEATURE,
            "droplet barcode": DROPLET_FEATURE, "UMI": UMI_FEATURE,
            "i7 sample index read": I7_FEATURE}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


T5_HEAD = il.TRUSEQ_READ2[14:]
ROUND1_BARCODE = "AAAGAA"
S5R_TAIL = "TCTTTCCCTACACGACGCTC"
BEAD_HEAD = "CTACACGACGCTCTTCCGATCT"
BEAD_LINKER = "TTTCTTATAT"
SAMPLE_INDEX_OLIGO = "CGTGAT"
SAMPLE_INDEX_READ = revcomp(SAMPLE_INDEX_OLIGO)


def indexed_tn5_top() -> list[Segment]:
    return [seg("TruSeq Read-2 head", T5_HEAD, "r2"),
            seg("round-1 barcode", "B" * 6, "cbc", placeholder=True),
            seg("mosaic end", nx.ME, "me")]


def bead_tso() -> list[Segment]:
    return [seg("TruSeq Read-1 head", BEAD_HEAD, "r1"),
            seg("droplet barcode", "C" * 16, "cbc", placeholder=True),
            seg("UMI", "N" * 10, "umi", placeholder=True),
            seg("TSO linker", BEAD_LINKER), seg("rGrGrG", "rGrGrG", placeholder=True)]


def hybrid_tagment_scene(insert_nt: int = 28) -> Scene:
    top = [seg("cDNA", "X" * insert_nt, placeholder=True),
           *indexed_tn5_top()]
    sc = Scene()
    sc.strand("cDNA", top, label="RNA/cDNA hybrid, cDNA strand")
    sc.anneal("ME bottom", [seg("ME reverse complement", nx.ME_RC, "me")],
              to="cDNA", pair=("ME reverse complement", "mosaic end"),
              label="5'-phosphorylated, 3'-ddC")
    sc.mark("cDNA", "mosaic end", "round-1 indexed Tn5 end")
    return sc


SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])
RUN_ROLES = ("Read 1", "Index 1 (i7)", "Read 2")


def final_library(insert_nt: int = 28) -> Construct:
    lib = Construct([
        seg("P5", il.P5, "p5"), seg("S5R-P5 tail", S5R_TAIL, "r1"),
        seg("remaining Read-1 site", BEAD_HEAD[-9:], "r1"),
        seg("droplet barcode", "C" * 16, "cbc", placeholder=True),
        seg("UMI", "N" * 10, "umi", placeholder=True),
        seg("TSO linker and GGG", BEAD_LINKER + "GGG"),
        seg("5-prime transcript", "X" * insert_nt, placeholder=True),
        seg("mosaic end reverse complement", nx.ME_RC, "me"),
        seg("round-1 barcode reverse complement", "B" * 6, "cbc", placeholder=True),
        seg("TruSeq Read-2 reverse complement", revcomp(il.TRUSEQ_READ2), "r2"),
        seg("i7 sample index read", SAMPLE_INDEX_READ, "cbc"),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="FIPRESCI sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS, required_roles=RUN_ROLES)
    if problems:
        raise ValueError("invalid final FIPRESCI construct: " + "; ".join(problems))
    return lib


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")


def protocol_graph() -> ProtocolGraph:
    source = "Defining source: https://doi.org/10.1186/s13059-023-02893-1"
    graph = ProtocolGraph("FIPRESCI", authority=source)
    cells = graph.input("cells", "fixed permeabilized cells containing RNA",
                        MaterialKind.BIOLOGICAL, primary=True)
    rt_primer = graph.input("rt_primer", "in-cell reverse-transcription primer",
                            MaterialKind.OLIGO, inferred=True,
                            source=source + " (primer sequence not represented on this page)")
    rt = graph.input("reverse_transcriptase", "reverse transcriptase",
                     MaterialKind.CATALYST)
    hybrids = graph.product("rna_cdna_hybrids", "RNA/cDNA hybrids with 3'-CCC cDNA tails",
                            MaterialKind.MOLECULAR)
    graph.react("reverse_transcription", "Reverse-transcribe RNA inside fixed cells",
                StepKind.REVERSE_TRANSCRIPTION,
                ((cells, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (hybrids,))
    indexed_adapter = graph.input("indexed_transposome", "round-1 indexed Tn5 adapter",
                                  MaterialKind.ADAPTER,
                                  state=tuple(indexed_tn5_top()))
    tn5 = graph.input("tn5", "Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagmented_cells", "cells containing indexed hybrid-tagmented cDNA",
                           MaterialKind.MOLECULAR, state=hybrid_tagment_scene())
    tagment_waste = graph.product("tagmentation_waste", "non-library tagmentation products",
                                  MaterialKind.WASTE)
    graph.react("hybrid_tagment", "Tagment RNA/cDNA hybrids in indexed wells",
                StepKind.TAGMENTATION,
                ((hybrids, InputRole.TEMPLATE), (indexed_adapter, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tagment_waste))
    graph.discard(tagment_waste,
                  "fragments lacking a usable indexed Tn5 end are not enriched")
    loading = graph.input("gem_loading", "overloaded 10x GEM partition condition",
                          MaterialKind.CONDITION)
    droplets = graph.product("productive_gems", "GEMs containing cells and gel-bead oligos",
                             MaterialKind.MOLECULAR)
    excluded = graph.product("excluded_gems", "empty or multiplet GEMs",
                             MaterialKind.WASTE)
    graph.react("partition", "Overload a 10x 5-prime channel with tagmented cells",
                StepKind.PARTITION,
                ((tagged, InputRole.TEMPLATE), (loading, InputRole.CONDITION)),
                (droplets, excluded))
    graph.discard(excluded, "empty and multiplet GEMs are excluded from accepted cells")
    bead = graph.input("bead_tso", "10x 5-prime gel-bead template-switch oligo",
                       MaterialKind.OLIGO, state=tuple(bead_tso()))
    switched = graph.product("dual_index_cdna", "plate- and droplet-indexed 5-prime cDNA",
                             MaterialKind.MOLECULAR)
    graph.react("template_switch", "Template-switch onto the released bead TSO",
                StepKind.TEMPLATE_SWITCHING,
                ((droplets, InputRole.TEMPLATE), (bead, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (switched,))
    forward = graph.input("enrichment_forward", "S5R-P5 enrichment primer",
                          MaterialKind.OLIGO)
    reverse = graph.input("enrichment_reverse", "Tn5-end enrichment primer",
                          MaterialKind.OLIGO)
    polymerase = graph.input("pcr_polymerase", "library PCR polymerase",
                             MaterialKind.CATALYST)
    enriched = graph.product("enriched_products", "TSO-bearing enriched cDNA products",
                             MaterialKind.MOLECULAR)
    graph.react("enrichment_pcr", "Select molecules bearing the droplet TSO and indexed Tn5 end",
                StepKind.PCR,
                ((switched, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (reverse, InputRole.STRUCTURAL), (polymerase, InputRole.CATALYST)),
                (enriched,), geometry=PCRGeometry.unavailable(
                    final_library(), reason="The pre-enrichment hybrid fragment is variable; the final primer sites are verified separately."))
    beads = graph.input("streptavidin_beads", "streptavidin capture beads",
                        MaterialKind.SELECTION)
    captured = graph.product("captured_products", "biotin-selected enrichment products",
                             MaterialKind.MOLECULAR)
    uncaptured = graph.product("uncaptured_products", "uncaptured amplification products",
                               MaterialKind.WASTE)
    graph.react("capture", "Capture the biotin-bearing enrichment products",
                StepKind.AFFINITY_CAPTURE,
                ((enriched, InputRole.TEMPLATE), (beads, InputRole.SELECTION)),
                (captured, uncaptured))
    graph.discard(uncaptured, "unbound products are removed during capture washes")
    sample_primer = graph.input("sample_index_primer", "S-P7 sample-index primer",
                                MaterialKind.OLIGO)
    library = graph.product("library", final_library().name, MaterialKind.LIBRARY,
                            state=final_library())
    graph.react("library_pcr", "Complete the adapters and add the sample index",
                StepKind.PCR,
                ((captured, InputRole.TEMPLATE), (forward, InputRole.STRUCTURAL),
                 (sample_primer, InputRole.STRUCTURAL),
                 (polymerase, InputRole.CATALYST)), (library,),
                geometry=PCRGeometry.unavailable(
                    final_library(), reason="The captured pre-library intermediate is not represented base-for-base; final sites are exact."))
    graph.sequence(library, SequencingPlan(
        library=final_library(), primers=SEQ_PRIMERS, acquired_roles=RUN_ROLES,
        read_lengths={"Read 1": 150, "Index 1 (i7)": 6, "Read 2": 150},
        read_uses={"Read 1": "extract_then_include",
                   "Index 1 (i7)": "demultiplex",
                   "Read 2": "extract_then_include"}))
    graph.declare_reviewed_scope((
        "cells", "rt_primer", "reverse_transcriptase", "reverse_transcription",
        "indexed_transposome", "tn5", "hybrid_tagment", "gem_loading",
        "partition", "bead_tso", "template_switch", "enrichment_forward",
        "enrichment_reverse", "pcr_polymerase", "enrichment_pcr",
        "streptavidin_beads", "capture", "sample_index_primer", "library_pcr",
        "library"), source)
    return graph.require_publishable()
