from types import SimpleNamespace
from protocol_model import chain, role_library
from multimodal_protocols import CELL, Definition, InputRole, MaterialKind, ProtocolGraph, SAMPLE, StepKind, _finish_pcr, ident, molecule, role, role_library
DEFINITION = Definition("CoTECH — chromatin occupancy and transcriptome", "10.1038/s41592-021-01129-z",
        "Xiong et al., Nature Methods (2021)", "Barcoded antibody-targeted pA–Tn5 and barcoded poly(dT) RT mark chromatin and RNA in the same first well; redistribution supplies a second combinatorial cell index.",
        "Separate chromatin and RNA libraries share the ordered pair of well barcodes.", "01_cotech.html")


def _cotech_graph(d, final):
    source = f"{d.citation}; doi:{d.doi}"
    graph = ProtocolGraph(d.title, authority=source)
    cells = graph.input("cells", "antibody-bound permeabilized cells containing poly(A) RNA",
                        MaterialKind.BIOLOGICAL, state=_cotech().initial, primary=True)
    transposome = graph.input("barcoded_patn5", "first-well-barcoded pA–Tn5",
                              MaterialKind.ADAPTER)
    tn5 = graph.input("tn5", "protein-A–Tn5 transposase", MaterialKind.CATALYST)
    tagged = graph.product("tagged_cells", "cells with target-proximal barcoded chromatin",
                           MaterialKind.MOLECULAR, state=cells.state)
    tn5_waste = graph.product("tn5_waste", "non-library tagmentation products", MaterialKind.WASTE)
    graph.react("tagment", "Tagment antibody-targeted chromatin",
                StepKind.TAGMENTATION,
                ((cells, InputRole.TEMPLATE), (transposome, InputRole.STRUCTURAL),
                 (tn5, InputRole.CATALYST)), (tagged, tn5_waste))
    graph.discard(tn5_waste, "does not retain a target-proximal barcoded fragment")
    rt_primer = graph.input("rt_primer", "matching first-well-barcoded poly(dT) primer",
                            MaterialKind.OLIGO)
    rt = graph.input("reverse_transcriptase", "reverse transcriptase", MaterialKind.CATALYST)
    paired = graph.product("paired_branches", "cells containing first-well-barcoded chromatin and cDNA",
                           MaterialKind.MOLECULAR,
                           state=molecule("paired intracellular branches", (
                               role("barcoded chromatin fragments"),
                               role("RT well barcode"), role("transcript cDNA"))))
    graph.react("reverse_transcribe", "Reverse-transcribe RNA with the matching well barcode",
                StepKind.REVERSE_TRANSCRIPTION,
                ((tagged, InputRole.TEMPLATE), (rt_primer, InputRole.STRUCTURAL),
                 (rt, InputRole.CATALYST)), (paired,))
    second_well = graph.input("second_well", "second-round combinatorial-index partition",
                              MaterialKind.CONDITION)
    indexed = graph.product("indexed_products", "two-well-indexed chromatin and cDNA products",
                            MaterialKind.MOLECULAR,
                            state=molecule("two-index cell identity", (
                                ident("first-round well barcode", "cotech_round1", "cell_barcode"),
                                ident("second-round well barcode", "cotech_round2", "cell_barcode"),
                                role("DNA and cDNA products"))))
    collisions = graph.product("collisions", "ambiguous second-well collisions", MaterialKind.WASTE)
    graph.react("redistribute", "Pool and redistribute products into second-index wells",
                StepKind.PARTITION,
                ((paired, InputRole.TEMPLATE), (second_well, InputRole.CONDITION)),
                (indexed, collisions), changes_molecule=False)
    graph.discard(collisions, "does not have an unambiguous ordered pair of cell indexes")
    chromatin = graph.product("chromatin_branch", "CoTECH chromatin branch",
                              MaterialKind.MOLECULAR, state=indexed.state)
    rna = graph.product("rna_branch", "CoTECH transcriptome branch",
                        MaterialKind.MOLECULAR, state=indexed.state)
    graph.react("split_libraries", "Split chromatin and transcriptome libraries",
                StepKind.SPLIT, ((indexed, InputRole.TEMPLATE),), (chromatin, rna),
                changes_molecule=False)
    reviewed = ["cells", "barcoded_patn5", "tn5", "tagment", "rt_primer",
                "reverse_transcriptase", "reverse_transcribe", "second_well",
                "redistribute", "split_libraries"]
    chrom_final = role_library("CoTECH chromatin-library roles",
                               (CELL, role("target-proximal chromatin insert")), index=SAMPLE)
    rna_final = role_library("CoTECH transcriptome-library roles",
                             (CELL, role("transcript cDNA")), index=SAMPLE)
    _finish_pcr(graph, "chromatin", chromatin, chrom_final,
                "Chromatin branch: " + d.endpoint, reviewed, label="CoTECH chromatin")
    _finish_pcr(graph, "rna", rna, rna_final,
                "RNA branch: " + d.endpoint, reviewed, label="CoTECH RNA")
    graph.declare_reviewed_scope(reviewed, source)
    return graph.require_publishable()


def protocol_graph():
    d = DEFINITION
    final = FINAL_LIBRARY
    primers = SEQ_PRIMERS
    key = 'cotech'
    return _cotech_graph(d, final)

def _cotech():
    return chain(molecule("antibody-bound permeabilized cells", (
        role("target chromatin epitope"), role("poly(A) RNA"))), (
        ("Bind uniquely barcoded protein-A–Tn5 and tagment target-proximal DNA",
         molecule("tagged chromatin branch", (role("T5 well barcode"),
                                               role("epitope-proximal DNA"),
                                               role("T7 well barcode"))),
         "The antibody localizes pA–Tn5; the first well identifies the chromatin branch."),
        ("Reverse-transcribe mRNA with a matching well-barcoded poly(dT) primer",
         molecule("paired intracellular branches", (role("barcoded chromatin fragments"),
                                                      role("RT well barcode"),
                                                      role("transcript cDNA"))),
         "Chromatin and RNA receive a shared first-round well identity before pooling."),
        ("Pool, redistribute and preamplify both branches with a second well index",
         molecule("two-index cell identity", (ident("first-round well barcode", "cotech_round1", "cell_barcode"),
                                               ident("second-round well barcode", "cotech_round2", "cell_barcode"),
                                               role("DNA and cDNA products"))),
         "The ordered pair of wells identifies a cell after combinatorial indexing."),
        ("Split the product into chromatin and transcriptome libraries",
         molecule("matched CoTECH libraries", (role("TruSeq-like chromatin library"),
                                                role("shared cell identity"),
                                                role("Nextera-like RNA library"))),
         "The two physical libraries remain linked by the same combinatorial barcode."),
    ))

d = DEFINITION
FINAL_LIBRARY = role_library(d.title + " final library roles", (CELL, role("chromatin or cDNA insert")), index=SAMPLE)
SPEC = SimpleNamespace(
    TITLE=d.title,
    NOTES=d.notes,
    SOURCE=f'Defining source: <a href="https://doi.org/{d.doi}">{d.citation}</a>.',
    SUMMARY=d.summary,
    CAVEAT="",
    FINAL_LIBRARY=FINAL_LIBRARY,
    SEQ_PRIMERS=(),
    FINAL_CAPTION=d.endpoint,
    SEQUENCING_INTRO=d.endpoint,
    READ_LENGTHS=None,
    SEQUENCING_ENDING=d.endpoint,
    SEQUENCING_UNAVAILABLE="Target panels, pool barcodes and kit revisions vary; the defining source establishes the read roles, while this page does not invent a universal primer sequence.",
    workflow=_cotech,
    oligos=lambda: (),
)
SPEC.protocol_graph = protocol_graph
globals().update(vars(SPEC))
