"""Molecular construct model for s3-WGS (Mulqueen et al., 2021)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Segment, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (Disposition, Terminal, InputRole,MaterialKind,ProtocolGraph,SequencingPlan,StepKind)


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


# Supplementary Tables 2--5.  One representative well/index combination is used in
# drawings; the variable regions are explicitly tagged as barcodes.
SBS12_PARTIAL = il.TRUSEQ_READ2[-18:]
TN5_INDEX = "GAACCGCG"                 # SBS12_18_UME_sci_1
I7_PRIMER_INDEX = "TCGCCTTA"           # PCR_i7_P7.S701, as carried in the oligo
I7_READ = revcomp(I7_PRIMER_INDEX)       # S701 / N701 bases reported by the sequencer
I5_INDEX = "CCTTAAGA"                   # PCR_A_i5_A, as carried in the oligo
DU = "U"

A14_LNA_BASES = nx.ADAPTOR_S5
A14_LNA_POSITIONS = (23, 25, 27, 29, 31, 33)  # one-based, as printed with + bases


def tn5_transfer_segments() -> list[Segment]:
    """Transferred strand of one indexed, U-containing Tn5 adaptor."""
    return [
        seg("partial TruSeq Read 2", SBS12_PARTIAL, "r2"),
        seg("Tn5 barcode", TN5_INDEX, "cbc"),
        seg("dU", DU, "w1", placeholder=True),
        seg("mosaic end", nx.ME, "me"),
    ]


def switching_oligo_segments() -> list[Segment]:
    return [seg("s5", nx.S5, "s5"), seg("LNA mosaic end", nx.ME, "me")]


def pre_switch_strand(insert_nt: int = 28) -> Construct:
    """One gap-filled strand: polymerase stopped immediately before dU."""
    return Construct([
        *tn5_transfer_segments(),
        seg("genomic insert", "X" * insert_nt, placeholder=True),
        seg("copied mosaic end", revcomp(nx.ME), "me"),
    ], name="gap-filled s3 strand")


def switched_strand(insert_nt: int = 28) -> Construct:
    """The same strand after A14-LNA-ME has templated addition of s5'."""
    return Construct([
        *pre_switch_strand(insert_nt),
        seg("copied s5", revcomp(nx.S5), "s5"),
    ], name="adapter-switched s3 strand")


def i7_pcr_primer() -> str:
    return il.P7 + I7_PRIMER_INDEX + il.TRUSEQ_READ2


def i5_pcr_primer() -> str:
    return nx.n5xx_primer(I5_INDEX)


def final_library(insert_nt: int = 28) -> Construct:
    """Final duplex represented by its P5-to-P7' strand (5' to 3')."""
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("i5", I5_INDEX, "cbc",
            feature=feature("cell_i5", "cell_barcode", "combinatorial", group="cell_id", part="i5")),
        seg("s5", nx.S5, "s5"),
        seg("mosaic end", nx.ME, "me"),
        seg("genomic insert", "X" * insert_nt, placeholder=True),
        seg("opposite mosaic end", nx.ME_RC, "me"),
        seg("copied dU", "A", "w1"),
        seg("Tn5 barcode, read orientation", revcomp(TN5_INDEX), "cbc",
            feature=feature("cell_tn5", "cell_barcode", "combinatorial", group="cell_id", part="Tn5")),
        seg("TruSeq Read 2 reverse complement", revcomp(il.TRUSEQ_READ2), "r2"),
        seg("i7", I7_READ, "cbc",
            feature=feature("cell_i7", "cell_barcode", "combinatorial", group="cell_id", part="i7")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="s3-WGS sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid final s3-WGS construct: " + "; ".join(problems))
    return lib


SEQ_PRIMERS = [
    sp.NEXTERA["R1"],
    sp.TRUSEQ["I1"],
    sp.NEXTERA["I2"],
    sp.TRUSEQ["R2"],
]
READ_LENGTHS = {"Read 1":85,"Index 1 (i7)":10,"Index 2 (i5)":10,"Read 2":85}


def protocol_graph() -> ProtocolGraph:
    source="Defining source: https://doi.org/10.1038/s41587-021-00962-z"
    graph=ProtocolGraph("s3-WGS",authority=source)
    genome=graph.input("nuclei","permeabilized single-cell genomic DNA",
                       MaterialKind.BIOLOGICAL,primary=True)
    transposome=graph.input("indexed_transposome","indexed U-ME transposome",
                            MaterialKind.ADAPTER,state=tuple(tn5_transfer_segments()))
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=graph.product("tagmented_fragments","symmetrically U-ME-tagmented genomic fragments",
                         MaterialKind.MOLECULAR)
    waste=graph.product("tagment_waste","candidate non-library tagmentation products",MaterialKind.MOLECULAR)
    graph.react("tagment","Tagment each first-round well with one indexed U-ME transposome",
                StepKind.TAGMENTATION,
                ((genome,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,waste))
    graph.mark_terminal(Terminal(waste.id,Disposition.UNRESOLVED,"Candidate insertion products and their physical recovery remain unresolved."))
    gap_pol=graph.input("gap_fill_polymerase","Nextera gap-fill polymerase",MaterialKind.CATALYST)
    dntps=graph.input("gap_fill_dntps","gap-fill dNTPs",MaterialKind.OLIGO)
    filled=graph.product("gap_filled",pre_switch_strand().name,MaterialKind.MOLECULAR,
                         state=pre_switch_strand())
    graph.react("gap_fill","Fill the Tn5 gap until polymerase stops at dU",StepKind.EXTENSION,
                ((tagged,InputRole.TEMPLATE),(dntps,InputRole.STRUCTURAL),
                 (gap_pol,InputRole.CATALYST)),(filled,))
    switch=graph.input("switching_oligo","A14-LNA-ME adapter-switching oligo",
                       MaterialKind.OLIGO,state=tuple(switching_oligo_segments()))
    switched=graph.product("switched",switched_strand().name,MaterialKind.MOLECULAR,
                           state=switched_strand())
    graph.react("adapter_switch","Template addition of s5: 10 adapter-switching cycles",
                StepKind.EXTENSION,
                ((filled,InputRole.TEMPLATE),(switch,InputRole.STRUCTURAL),
                 (gap_pol,InputRole.CATALYST)),(switched,))
    i5=graph.input("i5_primer","indexed i5 PCR primer",MaterialKind.OLIGO)
    i7=graph.input("i7_primer","indexed i7 PCR primer",MaterialKind.OLIGO)
    q5u=graph.input("q5u","uracil-tolerant Q5U polymerase",MaterialKind.CATALYST)
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final)
    graph.react("library_pcr","Complete flow-cell adapters: 13–15 Q5U PCR cycles",
                StepKind.PCR,
                ((switched,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
                 (i7,InputRole.STRUCTURAL),(q5u,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.exact(
                    switched_strand(),PCRPrimer("i7 PCR primer",i7_pcr_primer(),anneal_nt=len(SBS12_PARTIAL)),
                    PCRPrimer("i5 PCR primer",i5_pcr_primer(),anneal_nt=len(nx.S5)),
                    statement="Supplementary Tables 4-5 print the i5 and i7 primers and their s5/partial-Read-2 landing regions."))
    graph.sequence(library,SequencingPlan(
        library=final,primers=tuple(SEQ_PRIMERS),
        acquired_roles=("Read 1","Index 1 (i7)","Index 2 (i5)","Read 2"),
        read_lengths=READ_LENGTHS,
        read_uses={"Read 1":"include","Index 1 (i7)":"extract",
                   "Index 2 (i5)":"extract","Read 2":"extract_then_include"}))
    graph.declare_reviewed_scope((
        "nuclei","indexed_transposome","tn5","tagment","gap_fill_polymerase",
        "gap_fill_dntps","gap_fill","switching_oligo","adapter_switch","i5_primer",
        "i7_primer","q5u","library_pcr","library"),source)
    return graph.require_publishable()
