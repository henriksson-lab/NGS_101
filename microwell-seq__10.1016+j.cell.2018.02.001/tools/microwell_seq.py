"""Molecular construct model for Microwell-seq (Han et al., 2018)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import rt
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, revcomp


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


BEAD_HEAD = "TTTAGGGATAACAGGGTAAT"
BEAD_END = "ACGT"
LINK1 = "CGACTCACTACAGGG"
LINK2 = "TCGGTGACACGATCG"
P5_SPACER = "GCCTGTCCGCGG"
I7_OLIGO_INDEX = "TAAGGCGA"
I7_INDEX_READ = revcomp(I7_OLIGO_INDEX)
TSO = rt.SMART_HANDLE + "GAAT" + "GGG"
P5_PRIMER = il.P5 + P5_SPACER + rt.SMART_HANDLE + "AC"
CUSTOM_R1 = P5_SPACER + rt.SMART_HANDLE + BEAD_END


def barcode(name: str) -> Segment:
    return seg(name, "B" * 6, "cbc", placeholder=True)


def bead_oligo() -> Construct:
    """Published seqA bead architecture after three extension rounds."""
    return Construct([
        seg("bead head", BEAD_HEAD), seg("SMART handle", rt.SMART_HANDLE, "tso"),
        seg("seqA constant", BEAD_END), barcode("barcode 1"), seg("linker 1", LINK1),
        barcode("barcode 2"), seg("linker 2", LINK2), barcode("barcode 3"),
        seg("UMI", "N" * 6, "umi", placeholder=True), seg("poly(T)", "T" * 30),
    ], name="Microwell-seq seqA bead oligo")


def tso_oligo() -> list[Segment]:
    return [seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("GAAT", "GAAT"),
            seg("rGrG + modified G", "rGrG+G", placeholder=True)]


def rt_scene(insert_nt: int = 28) -> Scene:
    mrna = [seg("transcript", "X" * insert_nt, placeholder=True), seg("poly(A)", "A" * 30)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA", mod3="poly(A)")
    sc.anneal("bead", list(bead_oligo()), to="mRNA", pair=("poly(T)", "poly(A)"),
              label="barcoded bead primer")
    sc.arrow("bead", "reverse transcription and template switching")
    return sc


def inferred_s7_end() -> list[Segment]:
    """Same-lab/primer-supported, but not printed for the defining 2018 protocol."""
    return [seg("mosaic end reverse complement", nx.ME_RC, "me", inferred=True),
            seg("s7 reverse complement", nx.S7_RC, "s7", inferred=True),
            seg("i7 index read", I7_INDEX_READ, "cbc", inferred=True),
            seg("P7 reverse complement", il.P7_RC, "p7", inferred=True)]


SEQ_PRIMERS = (
    sp.custom("Read 1", "Microwell-seq custom Read 1", CUSTOM_R1,
              "Han et al. 2018, Supplementary Table S2"),
    sp.NEXTERA["I1"], sp.NEXTERA["R2"],
)
RUN_ROLES = ("Read 1", "Index 1 (i7)", "Read 2")


def final_library(insert_nt: int = 28) -> Construct:
    bead = list(bead_oligo())
    lib = Construct([
        seg("P5", il.P5, "p5"), seg("custom spacer", P5_SPACER),
        *bead[1:], seg("cDNA", "X" * insert_nt, placeholder=True), *inferred_s7_end(),
    ], name="Microwell-seq sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS, required_roles=RUN_ROLES)
    if problems:
        raise ValueError("invalid final Microwell-seq construct: " + "; ".join(problems))
    return lib


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")


def read1_layout() -> list[tuple[str, str]]:
    return [("1&ndash;6", "cell barcode 1"), ("7&ndash;21", "linker 1"),
            ("22&ndash;27", "cell barcode 2"), ("28&ndash;42", "linker 2"),
            ("43&ndash;48", "cell barcode 3"), ("49&ndash;54", "UMI")]
