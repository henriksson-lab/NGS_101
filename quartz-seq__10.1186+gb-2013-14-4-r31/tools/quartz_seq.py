"""Molecular construct model for Quartz-Seq (Sasagawa et al., 2013)."""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, revcomp


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


M = "TATAGAATTCGCGGCCGCTCGCGAT"
RT_MIDDLE = "AATACGACTCACTATAGGGCG"
RT_PRIMER = M + RT_MIDDLE + "T" * 24
TAGGING_PRIMER = M + "T" * 23
SUPPRESSION_PRIMER = "G" + M
TRSU = il.TRUSEQ_P5_FULL
TRSI2_INDEX = "CGATGT"
TRSI2 = il.INDEX1_PRIMER + TRSI2_INDEX + il.P7_RC
TPC1 = il.P5[:21]
TPC2 = il.P7[:22]


def rt_primer() -> list[Segment]:
    return [seg("M handle", M, "tso"), seg("T7/GGCG", RT_MIDDLE),
            seg("dT24", "T" * 24)]


def tagging_primer() -> list[Segment]:
    return [seg("M handle", M, "tso"), seg("dT23", "T" * 23)]


def mrna_rt_scene() -> Scene:
    mrna = [seg("RNA body", "X" * 28, placeholder=True), seg("poly(A)", "A" * 24)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA")
    sc.anneal("RT primer", rt_primer(), to="mRNA", pair=("dT24", "poly(A)"),
              label="RT primer")
    sc.arrow("RT primer", "reverse transcriptase")
    return sc


def tailed_first_strand() -> Construct:
    return Construct([
        *rt_primer(),
        seg("cDNA", "X" * 28, placeholder=True),
        seg("TdT poly(A)", "A" * 24),
    ], name="poly(A)-tailed Quartz-Seq first strand")


def tagging_scene() -> Scene:
    target = tailed_first_strand()
    sc = Scene()
    sc.strand("first strand", list(target), label="first strand")
    sc.anneal("tagging primer", tagging_primer(), to="first strand",
              pair=("dT23", "TdT poly(A)"), label="tagging primer")
    sc.arrow("tagging primer", "second-strand synthesis")
    return sc


def wta_amplicon() -> Construct:
    return Construct([
        seg("tagging-end M", M, "tso"),
        seg("tagging poly(T)", "T" * 23),
        seg("cDNA", "X" * 32, placeholder=True),
        seg("RT-primer poly(T)", "T" * 24),
        seg("T7/GGCG reverse complement", revcomp(RT_MIDDLE)),
        seg("RT-end M reverse complement", revcomp(M), "tso"),
    ], name="Quartz-Seq WTA product")


SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])


def final_library(insert_nt: int = 36) -> Construct:
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("TruSeq Read 1 remainder", il.TRUSEQ_READ1[4:], "r1"),
        seg("cDNA fragment", "X" * insert_nt, placeholder=True),
        seg("TruSeq Read 2 reverse complement", revcomp(il.TRUSEQ_READ2), "r2"),
        seg("i7", TRSI2_INDEX, "cbc"),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="Quartz-Seq sequencing library")
    for primer in SEQ_PRIMERS:
        if sp.locate(lib, primer) is None:
            raise ValueError(f"{primer.name} does not land on {lib.name}")
    return lib


def primer_landings(lib: Construct | None = None):
    lib = lib or final_library()
    out = []
    for primer in SEQ_PRIMERS:
        hit = sp.locate(lib, primer)
        if hit is None:
            raise ValueError(f"{primer.name} does not land on {lib.name}")
        out.append((primer, hit))
    return out
