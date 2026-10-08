"""Molecular model for the original Drop-seq implementation of CROP-seq."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import nextera as nx
import seqprimers as sp
from crispr import SCAFFOLD_V1
from chemdraw import Construct, Segment, revcomp
from lentivirus import TransferGenome


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def role(name: str, text: str, tag: str | None = None) -> Segment:
    return seg(name, f"[{text}]", tag, placeholder=True)


GIBSON_5 = "TGGAAAGGACGAAACACCG"
GIBSON_3 = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGC"
GUIDE_LEN = 20
SMART_HANDLE = "AAGCAGTGGTATCAACGCAGAGT"
P5_SPACER = "GCCTGTCCGCGG"
CUSTOM_R1 = P5_SPACER + SMART_HANDLE + "AC"


def pooled_guide_oligo() -> list[Segment]:
    return [seg("hU6 homology + G", GIBSON_5, "r3"),
            seg("guide", "N" * GUIDE_LEN, "cbc", placeholder=True),
            seg("scaffold homology", GIBSON_3, "r2")]


def guide_expression_cassette() -> tuple[Segment, ...]:
    return (seg("hU6 promoter tail + G", GIBSON_5, "r3"),
            seg("guide", "N" * GUIDE_LEN, "cbc", placeholder=True),
            seg("gRNA scaffold", SCAFFOLD_V1, "r2"))


def transfer_genome() -> TransferGenome:
    return TransferGenome(
        u3=(*guide_expression_cassette(), role("U3 remainder", "U3 remainder", "me")),
        r=(role("R", "R", "me"),), u5=(role("U5", "U5", "me"),),
        internal=(role("vector body", "psi ... EF1a-Puro ... WPRE", "w1"),),
        name="CROP-seq")


def pol3_guide_transcript() -> Construct:
    return Construct([seg("guide", "N" * GUIDE_LEN, "cbc", placeholder=True),
                      seg("gRNA scaffold", SCAFFOLD_V1, "r2"),
                      role("Pol III termination", "U6 termination")],
                     name="functional Pol III guide RNA")


def pol2_reporter_transcript() -> Construct:
    return Construct([role("puromycin transcript", "EF1a--Puro--WPRE", "w1"),
                      *guide_expression_cassette(), role("3' LTR", "3' LTR", "me"),
                      seg("poly(A)", "A" * 24)],
                     name="polyadenylated Pol II reporter RNA")


def guide_bearing_cdna_fragment() -> Construct:
    """Representative bead-end fragment cut at the hU6/guide boundary."""
    return Construct([
        seg("SMART handle", SMART_HANDLE, "tso"), seg("batch-B constant", "AC"),
        seg("cell barcode", "B" * 12, "cbc", placeholder=True),
        seg("UMI", "U" * 8, "umi", placeholder=True), seg("poly(T)", "T" * 30),
        seg("3' transcript cDNA", "X" * 16, "me", placeholder=True),
        seg("scaffold cDNA", revcomp(SCAFFOLD_V1), "r2"),
        seg("guide cDNA", "N" * GUIDE_LEN, "cbc", placeholder=True),
        seg("mosaic end reverse complement", nx.ME_RC, "me"),
        seg("s7 reverse complement", nx.S7_RC, "s7")],
        name="guide-bearing bead-end fragment")


def sequencing_primers():
    return (sp.custom("Read 1", "Drop-seq Custom Read1", CUSTOM_R1,
                      "CROP-seq Supplementary Protocol oligo-order table"),
            sp.NEXTERA["I1"], sp.NEXTERA["R2"])


def final_library() -> Construct:
    lib = Construct([seg("P5", il.P5, "p5"), seg("Read-1 spacer", P5_SPACER, "r1"),
                     *list(guide_bearing_cdna_fragment()),
                     seg("i7 reverse complement", "I" * 8, "cbc", placeholder=True),
                     seg("P7 reverse complement", il.P7_RC, "p7")],
                    name="CROP-seq guide-bearing Drop-seq library")
    problems = sp.verify(lib, sequencing_primers(),
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if problems:
        raise ValueError("invalid CROP-seq final library: " + "; ".join(problems))
    return lib


def read1_layout() -> list[tuple[str, str]]:
    return [("1-12", "cell barcode"), ("13-20", "UMI")]
