"""Vitak et al. 2017 SCI-seq four-index library geometry."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Row, Scene, feature, revcomp
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

TITLE = "SCI-seq — single-cell combinatorial-indexed genome sequencing"
NOTES = "01_sci-seq.html"
SOURCE = ('Defining paper: <a href="https://doi.org/10.1038/nmeth.4154">Vitak et al., '
          '<i>Nature Methods</i> (2017)</a>; inherited indexed-transposome oligos: '
          '<a href="https://doi.org/10.1038/ng.3119">Amini et al. (2014)</a>.')
SUMMARY = ("Nucleosome-depleted nuclei receive a first bipartite barcode during indexed "
           "transposition, are pooled and resorted, then receive a second bipartite barcode "
           "during PCR. The four index components identify one nucleus.")
CAVEAT = "Representative index bases are placeholders; all constant connector and primer bases are source-transcribed."

# Amini Supplementary Table 4; SCI-seq explicitly reuses these complexes and primers.
C15 = "GCGATCGAGGACGGC"
D15 = "CACCGTCTCCGCCTC"
CONNECTOR_A = "TCCACGC"
CONNECTOR_B = "CTGTCCCTGTCC"
R1_SEQ = C15 + nx.ME
R2_SEQ = D15 + nx.ME
I1_SEQ = nx.ME_RC + revcomp(D15)
I2_SEQ = nx.ME_RC + revcomp(C15)

FINAL_LIBRARY = Construct([
    seg("P5", il.P5, "p5"),
    seg("PCR i5", "J" * 10, "cbc", placeholder=True,
        feature=feature("cell_pcr_i5", "cell_barcode", "combinatorial", group="cell_id", part="PCR i5")),
    seg("S5", nx.S5, "s5"), seg("A connector", CONNECTOR_A, "me"),
    seg("transposase i5", "K" * 8, "cbc", placeholder=True,
        feature=feature("cell_tn5_i5", "cell_barcode", "combinatorial", group="cell_id", part="transposase i5")),
    seg("C15", C15, "me"), seg("left mosaic end", nx.ME, "me"),
    seg("genomic DNA", "X" * 38, placeholder=True),
    seg("right mosaic end reverse complement", nx.ME_RC, "me"),
    seg("D15 reverse complement", revcomp(D15), "me"),
    seg("transposase i7 reverse complement", "L" * 8, "cbc", placeholder=True,
        feature=feature("cell_tn5_i7", "cell_barcode", "combinatorial", group="cell_id", part="transposase i7")),
    seg("B connector reverse complement", revcomp(CONNECTOR_B), "me"),
    seg("S7 reverse complement", nx.S7_RC, "s7"),
    seg("PCR i7 reverse complement", "I" * 10, "cbc", placeholder=True,
        feature=feature("cell_pcr_i7", "cell_barcode", "combinatorial", group="cell_id", part="PCR i7")),
    seg("P7 reverse complement", il.P7_RC, "p7"),
], name="quad-index SCI-seq library")

SEQ_PRIMERS = (
    sp.custom("Read 1", "SCI-seq Read 1", R1_SEQ, "Amini et al. 2014 Supplementary Table 4"),
    sp.custom("Index 1 (i7)", "SCI-seq Index 1", I1_SEQ, "Amini et al. 2014 Supplementary Table 4"),
    sp.custom("Index 2 (i5)", "SCI-seq Index 2", I2_SEQ,
              "reverse-complement counterpart of the source-defined C15/ME Read 1 site"),
    sp.custom("Read 2", "SCI-seq Read 2", R2_SEQ, "Amini et al. 2014 Supplementary Table 4"),
)
if problems := sp.verify(FINAL_LIBRARY, SEQ_PRIMERS):
    raise ValueError("SCI-seq final library: " + "; ".join(problems))

FINAL_CAPTION = ("One productive heterologous molecule. Each index read acquires an 8-base "
                 "transposase index, skips constant connector bases, then acquires a 10-base PCR index.")
SEQUENCING_INTRO = ("The custom recipe is R1 50 cycles; I1 8 imaged, 27 dark, 10 imaged; "
                    "I2 8 imaged, 21 dark, 10 imaged; R2 50 cycles.")


def tagged_scene() -> Scene:
    con = Construct(FINAL_LIBRARY.segments[3:12], name="indexed transposition product")
    sc = Scene.duplex(list(con), label="inside one intact nucleus")
    sc.junction("top", "left mosaic end", "genomic DNA", "Tn5 transfer")
    sc.junction("top", "genomic DNA", "right mosaic end reverse complement", "Tn5 transfer")
    sc.labels("top")
    return sc


def sections():
    return [
        ("Expose genomic DNA while retaining nuclei", [
            Row(chunks=[("intact nucleus: histone-bound DNA -> LAND or crosslink/SDS nucleosome depletion", None, False)]),
        ], "Removing nucleosomes makes transposition genome-wide rather than accessibility-biased."),
        ("Indexed transposition in 96 wells", tagged_scene().rows(),
         "Each well supplies one i5/i7 transposome combination; ** marks the two Tn5-transfer boundaries."),
        ("Pool and resort nuclei", [
            Row(chunks=[("96 transposition wells -> pool -> FANS: 22 nuclei into each PCR well", None, False)]),
        ], "The nucleus keeps its first index pair while redistribution creates the second combinatorial dimension."),
        ("Add the PCR-well index pair", [
            Row(chunks=[("72 °C gap extension -> real-time PCR with one 10-base i5 and one 10-base i7", "cbc", False)]),
        ], "One cell is identified by transposase i5 + transposase i7 + PCR i5 + PCR i7."),
    ]


def protocol_graph():
    graph=ProtocolGraph(TITLE,authority=SOURCE)
    nuclei=graph.input("nuclei","intact nuclei",MaterialKind.BIOLOGICAL,primary=True)
    depletion=graph.input("depletion_condition","LAND or crosslink/SDS nucleosome depletion",
                           MaterialKind.CONDITION)
    depleted=graph.product("depleted","nucleosome-depleted intact nuclei",
                           MaterialKind.BIOLOGICAL)
    graph.react("deplete_nucleosomes","Expose genomic DNA while retaining nuclei",
                StepKind.GENERIC,
                ((nuclei,InputRole.TEMPLATE),(depletion,InputRole.CONDITION)),
                (depleted,),changes_molecule=False)
    transposome=graph.input("indexed_transposome","well-indexed i5/i7 Tn5 transposomes",
                            MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=graph.product("tagged","indexed transposition products in intact nuclei",
                         MaterialKind.MOLECULAR,state=tagged_scene())
    tagment_waste=graph.product("tagment_waste","non-amplifiable transposition products",
                                MaterialKind.WASTE)
    graph.react("tagment","Install the first bipartite index pair by transposition",
                StepKind.TAGMENTATION,
                ((depleted,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,tagment_waste))
    graph.discard(tagment_waste,"non-amplifiable Tn5 outcomes do not enter the library")
    sorting=graph.input("sorting_rule","FANS redistribution of 22 nuclei per PCR well",
                        MaterialKind.SELECTION)
    resorted=graph.product("resorted","pooled and resorted indexed nuclei",
                           MaterialKind.MOLECULAR)
    excluded=graph.product("excluded","nuclei excluded during sorting",MaterialKind.WASTE)
    graph.react("resort","Pool and resort nuclei into PCR wells",StepKind.SAMPLING,
                ((tagged,InputRole.TEMPLATE),(sorting,InputRole.SELECTION)),
                (resorted,excluded),changes_molecule=False)
    graph.discard(excluded,"nuclei outside the sorted PCR-well allocation are not processed")
    fill=graph.input("gap_polymerase","Tn5 gap-extension polymerase",MaterialKind.CATALYST)
    extended=graph.product("extended","gap-extended indexed fragments",MaterialKind.MOLECULAR)
    graph.react("gap_extend","Extend the transposase-created gaps",StepKind.EXTENSION,
                ((resorted,InputRole.TEMPLATE),(fill,InputRole.CATALYST)),(extended,))
    i5=graph.input("pcr_i5","10-base PCR-well i5 primer",MaterialKind.OLIGO)
    i7=graph.input("pcr_i7","10-base PCR-well i7 primer",MaterialKind.OLIGO)
    polymerase=graph.input("pcr_polymerase","real-time library PCR polymerase",
                           MaterialKind.CATALYST)
    library=graph.product("library",FINAL_LIBRARY.name,MaterialKind.LIBRARY,
                          state=FINAL_LIBRARY)
    graph.react("library_pcr","Add the second bipartite PCR-well index pair",
                StepKind.PCR,
                ((extended,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
                 (i7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),
                (library,),geometry=PCRGeometry.unavailable(
                    FINAL_LIBRARY,reason="The final four-index geometry is exact, while individual well index members vary."))
    graph.sequence(library,SequencingPlan(library=FINAL_LIBRARY,
                                          primers=tuple(SEQ_PRIMERS)))
    graph.declare_reviewed_scope((
        "nuclei","depletion_condition","deplete_nucleosomes","indexed_transposome",
        "tn5","tagment","sorting_rule","resort","gap_polymerase","gap_extend",
        "pcr_i5","pcr_i7","pcr_polymerase","library_pcr","library"),SOURCE)
    return graph.require_publishable()
