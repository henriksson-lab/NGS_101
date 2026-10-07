"""Molecular construct model for s3-WGS (Mulqueen et al., 2021)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Segment, revcomp


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
        seg("i5", I5_INDEX, "cbc"),
        seg("s5", nx.S5, "s5"),
        seg("mosaic end", nx.ME, "me"),
        seg("genomic insert", "X" * insert_nt, placeholder=True),
        seg("opposite mosaic end", nx.ME_RC, "me"),
        seg("copied dU", "A", "w1"),
        seg("Tn5 barcode, read orientation", revcomp(TN5_INDEX), "cbc"),
        seg("TruSeq Read 2 reverse complement", revcomp(il.TRUSEQ_READ2), "r2"),
        seg("i7", I7_READ, "cbc"),
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
