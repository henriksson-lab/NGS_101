"""Molecular construct model for Seq-Well S3 (Hughes et al., 2020)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import rt
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, revcomp


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


BEAD_PREFIX = "T" * 7
BEAD_SUFFIX = "AC"
P5_SPACER = "GCCTGTCCGCGG"
S3_RANDOM_TAIL = "NNNGGNNNB"
I7_INDEX_READ = "TAAGGCGA"
TSO = rt.SMART_HANDLE + "GAAT" + "GGG"
P5_HYBRID = il.P5 + P5_SPACER + rt.SMART_HANDLE + BEAD_SUFFIX
CUSTOM_R1 = P5_SPACER + rt.SMART_HANDLE + BEAD_SUFFIX


def bead_oligo() -> Construct:
    return Construct([
        seg("bead prefix", BEAD_PREFIX), seg("SMART handle", rt.SMART_HANDLE, "tso"),
        seg("bead constant", BEAD_SUFFIX),
        seg("cell barcode", "B" * 12, "cbc", placeholder=True),
        seg("UMI", "N" * 8, "umi", placeholder=True), seg("poly(T)", "T" * 30),
    ], name="Seq-Well S3 bead oligo")


def tso_oligo() -> list[Segment]:
    return [seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("GAAT", "GAAT"),
            seg("rGrGrG", "rGrGrG", placeholder=True)]


def randomer() -> list[Segment]:
    return [seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("GA", "GA"),
            seg("random priming tail", S3_RANDOM_TAIL, placeholder=True)]


def rt_scene(insert_nt: int = 28) -> Scene:
    mrna = [seg("transcript", "X" * insert_nt, placeholder=True),
            seg("poly(A)", "A" * 30)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA", mod3="poly(A)")
    sc.anneal("bead", list(bead_oligo()), to="mRNA", pair=("poly(T)", "poly(A)"),
              label="bead primer")
    sc.arrow("bead", "reverse transcription; template switching may occur")
    return sc


def random_second_strand_scene() -> Scene:
    # Placeholder text is deliberately symmetric under the Scene's reverse placement:
    # the exact random landing sequence is molecule-specific and is not asserted here.
    template = [seg("first-strand cDNA", "X" * 20, placeholder=True),
                seg("randomer landing", "BNNNCCNNN", placeholder=True)]
    sc = Scene()
    sc.strand("first", template, label="bead-bound first strand")
    sc.anneal("randomer", randomer(), to="first",
              pair=("random priming tail", "randomer landing"), label="S3 randomer")
    sc.arrow("randomer", "Klenow exo-minus makes the second strand toward the bead")
    return sc


SEQ_PRIMERS = (
    sp.custom("Read 1", "Seq-Well S3 custom Read 1", CUSTOM_R1,
              "Hughes et al. 2020 Key Resources Table"),
    sp.NEXTERA["I1"], sp.NEXTERA["R2"],
)
RUN_ROLES = ("Read 1", "Index 1 (i7)", "Read 2")


def final_library(insert_nt: int = 28) -> Construct:
    lib = Construct([
        seg("P5", il.P5, "p5"), seg("custom spacer", P5_SPACER),
        seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("bead constant", BEAD_SUFFIX),
        seg("cell barcode", "B" * 12, "cbc", placeholder=True),
        seg("UMI", "N" * 8, "umi", placeholder=True), seg("poly(T)", "T" * 30),
        seg("cDNA", "X" * insert_nt, placeholder=True),
        seg("mosaic end reverse complement", nx.ME_RC, "me"),
        seg("s7 reverse complement", nx.S7_RC, "s7"),
        seg("i7 index read", I7_INDEX_READ, "cbc"),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="Seq-Well S3 sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS, required_roles=RUN_ROLES)
    if problems:
        raise ValueError("invalid final Seq-Well S3 construct: " + "; ".join(problems))
    return lib


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")
