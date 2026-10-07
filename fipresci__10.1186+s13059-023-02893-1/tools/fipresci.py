"""Molecular construct model for FIPRESCI (Luo et al., 2023)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, revcomp


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
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
