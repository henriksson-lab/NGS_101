"""Molecular construct model for scifi-ATAC-seq (Wang et al.)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, revcomp


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


SPACER_A = "GATATGTGATAATGAGGAC"
SPACER_B = "TGAATGTGTAGAAGACAGA"
BC_A_EXAMPLE = "CGTAT"
BC_B_EXAMPLE = "TCGGA"
I7_INDEX_READ = "TAAGGCGA"
R1_SEQ = nx.S5 + SPACER_A
R2_SEQ = nx.S7 + SPACER_B
I1_SEQ = revcomp(R2_SEQ)
I2_SEQ = revcomp(R1_SEQ)


def adapter_a() -> list[Segment]:
    return [seg("s5", nx.S5, "s5"), seg("spacer A", SPACER_A),
            seg("Tn5 barcode A", "A" * 5, "cbc", placeholder=True),
            seg("mosaic end", nx.ME, "me")]


def adapter_b() -> list[Segment]:
    return [seg("s7", nx.S7, "s7"), seg("spacer B", SPACER_B),
            seg("Tn5 barcode B", "B" * 5, "cbc", placeholder=True),
            seg("mosaic end", nx.ME, "me")]


def bead_oligo() -> list[Segment]:
    """Unmodified 10x scATAC v1.1 gel-bead oligo from its vendor guide."""
    return [seg("P5", il.P5, "p5"),
            seg("GEM barcode", "C" * 16, "cbc", placeholder=True),
            seg("s5", nx.S5, "s5")]


def indexed_fragment(insert_nt: int = 28) -> Construct:
    return Construct([*adapter_a(), seg("accessible DNA", "X" * insert_nt, placeholder=True),
                      seg("opposite ME reverse complement", nx.ME_RC, "me"),
                      seg("Tn5 barcode B reverse complement", "B" * 5, "cbc", placeholder=True),
                      seg("spacer B reverse complement", revcomp(SPACER_B)),
                      seg("s7 reverse complement", nx.S7_RC, "s7")],
                     name="A/B-indexed accessible fragment")


SEQ_PRIMERS = (
    sp.custom("Read 1", "scifi 1_Read1", R1_SEQ, "Wang et al. Table S1"),
    sp.custom("Index 1 (i7)", "scifi 2_Index1(i7)", I1_SEQ, "Wang et al. Table S1"),
    sp.custom("Index 2 (i5)", "scifi 3_Index2(i5)", I2_SEQ, "Wang et al. Table S1"),
    sp.custom("Read 2", "scifi 4_Read2", R2_SEQ, "Wang et al. Table S1"),
)


def final_library(insert_nt: int = 28) -> Construct:
    lib = Construct([
        *bead_oligo(), seg("spacer A", SPACER_A),
        seg("Tn5 barcode A", "A" * 5, "cbc", placeholder=True),
        seg("mosaic end", nx.ME, "me"),
        seg("accessible genomic DNA", "X" * insert_nt, placeholder=True),
        seg("opposite ME reverse complement", nx.ME_RC, "me"),
        seg("Tn5 barcode B reverse complement", "B" * 5, "cbc", placeholder=True),
        seg("spacer B reverse complement", revcomp(SPACER_B)),
        seg("s7 reverse complement", nx.S7_RC, "s7"),
        seg("i7 sample index read", I7_INDEX_READ, "cbc"),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="scifi-ATAC-seq sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid final scifi-ATAC-seq construct: " + "; ".join(problems))
    return lib


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")
