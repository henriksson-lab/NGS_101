"""Molecular construct model for HyDrop-RNA (De Rop et al., 2022)."""
from __future__ import annotations

import illumina as il
import rt
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, revcomp


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


T7_PROMOTER = "TAATACGACTCACTATAGGG"
ACRYDITE_PRIMER = "T" * 7 + T7_PROMOTER + rt.SMART_HANDLE + "AC"
LINK1 = "CAGCTACTGC"
LINK2 = "CGAGTACCCT"
TSO = rt.SMART_HANDLE + "GAAT" + "GGG"
TSO_P = rt.SMART_HANDLE
I7_OLIGO_INDEX = "CGCTCAGTTC"
I5_OLIGO_INDEX = "TCGTGGAGCG"
HYI7_SPACER = "CTGTCCGCGG"
HYI7_SITE = HYI7_SPACER + rt.SMART_HANDLE + "AC"
CUSTOM_INDEX1 = revcomp(HYI7_SITE)
CUSTOM_READ2 = HYI7_SITE


def barcode(name: str) -> Segment:
    return seg(name, "B" * 10, "cbc", placeholder=True)


def bead_oligo() -> Construct:
    """Final RNA bead oligo; linkage and barcode order are fixed by construction."""
    return Construct([
        seg("poly(T) before promoter", "T" * 7),
        seg("T7 promoter", T7_PROMOTER),
        seg("SMART handle", rt.SMART_HANDLE, "tso"),
        seg("AC", "AC"),
        barcode("barcode 1"), seg("linker 1", LINK1),
        barcode("barcode 2"), seg("linker 2", LINK2),
        barcode("barcode 3"),
        seg("UMI", "N" * 8, "umi", placeholder=True),
        seg("poly(T)", "T" * 30),
    ], name="HyDrop-RNA bead oligo")


def tso_oligo() -> list[Segment]:
    return [seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("GAAT", "GAAT"),
            seg("rGrGrG", "rGrGrG", placeholder=True)]


def hyi7_primer() -> list[Segment]:
    return [seg("P7", il.P7, "p7"), seg("i7", I7_OLIGO_INDEX, "cbc"),
            seg("spacer", HYI7_SPACER), seg("SMART handle + AC", rt.SMART_HANDLE + "AC", "tso")]


def hyi5_primer() -> list[Segment]:
    return [seg("P5", il.P5, "p5"), seg("i5", I5_OLIGO_INDEX, "cbc"),
            seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1")]


def rt_scene(insert_nt: int = 28) -> Scene:
    """Released bead primer annealed to a poly(A)-tailed transcript."""
    primer = list(bead_oligo())
    mrna = [seg("transcript", "X" * insert_nt, placeholder=True),
            seg("poly(A)", "A" * 30)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA", mod3="poly(A)")
    sc.anneal("bead primer", primer, to="mRNA", pair=("poly(T)", "poly(A)"),
              label="released bead primer")
    sc.arrow("bead primer", "reverse transcription and template switching")
    return sc


def amplified_cdna(insert_nt: int = 28) -> Construct:
    """TSO-P amplicon; upstream acrydite/T7 prefix is excluded by the priming site."""
    bead = bead_oligo()
    keep = list(bead)[2:]
    return Construct([
        *keep,
        seg("cDNA", "X" * insert_nt, placeholder=True),
        seg("template-switch junction", "CCC"),
        seg("copied TSO suffix", revcomp("GAAT")),
        seg("opposite SMART handle", revcomp(rt.SMART_HANDLE), "tso"),
    ], name="HyDrop-RNA amplified cDNA")


def final_library(insert_nt: int = 28) -> Construct:
    """P5-to-P7' bead-end library selected by HYi5/HYi7 PCR."""
    lib = Construct([
        seg("P5", il.P5, "p5"), seg("i5 sample index", I5_OLIGO_INDEX, "cbc"),
        seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
        seg("cDNA", "X" * insert_nt, placeholder=True), seg("poly(A)", "A" * 30),
        seg("UMI", "N" * 8, "umi", placeholder=True), barcode("barcode 3"),
        seg("linker 2, plate orientation", revcomp(LINK2)), barcode("barcode 2"),
        seg("linker 1, plate orientation", revcomp(LINK1)), barcode("barcode 1"),
        seg("HYi7 site reverse complement", revcomp(HYI7_SITE), "tso"),
        seg("i7 reverse complement", revcomp(I7_OLIGO_INDEX), "cbc"),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="HyDrop-RNA sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid final HyDrop-RNA construct: " + "; ".join(problems))
    return lib


SEQ_PRIMERS = (
    sp.TRUSEQ["R1"],
    sp.custom("Index 1 (i7)", "HyDrop custom Index 1", CUSTOM_INDEX1,
              "De Rop et al. 2022, Supplementary file 3"),
    sp.TRUSEQ["I2"],
    sp.custom("Read 2", "HyDrop custom Read 2", CUSTOM_READ2,
              "De Rop et al. 2022, Supplementary file 3"),
)


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")


def read2_layout() -> list[tuple[str, str]]:
    names = [("barcode 1", "cell barcode 1"),
             ("linker 1, plate orientation", "linker 1"),
             ("barcode 2", "cell barcode 2"),
             ("linker 2, plate orientation", "linker 2"),
             ("barcode 3", "cell barcode 3"), ("UMI", "UMI")]
    rows, cycle = [], 1
    lib = final_library()
    for name, label in names:
        end = cycle + len(lib.get(name)) - 1
        rows.append((f"{cycle}&ndash;{end}", label))
        cycle = end + 1
    return rows
