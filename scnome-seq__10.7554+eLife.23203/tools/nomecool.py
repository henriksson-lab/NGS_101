"""Constructors for scNOMe-seq and scCOOL-seq library molecules."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Segment, feature, revcomp
import illumina as il
import seqprimers as sp

SCNOME = "scNOMe-seq"
SCCOOL = "scCOOL-seq"

RANDOM_NT = 9
COOL_P1_HANDLE = il.TRUSEQ_READ1[-22:]
COOL_P2_HANDLE = il.TRUSEQ_READ2[-21:]
FORWARD_P5 = il.P5 + il.TRUSEQ_READ1[4:]
INDEX_NT_NOME = 6
INDEX_NT_COOL = 8       # standard indexed-primer-shaped model; exact NEB set not printed

SOURCE_FORWARD_P5 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"
SOURCE_COOL_P1 = "CTACACGACGCTCTTCCGATCTNNNNNNNNN"


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def primer1_segments() -> list[Segment]:
    return [seg("Read 1 handle", COOL_P1_HANDLE, "r1"),
            seg("N9", "N" * RANDOM_NT, placeholder=True)]


def primer2_segments() -> list[Segment]:
    return [seg("Read 2 handle", COOL_P2_HANDLE, "r2"),
            seg("N9", "N" * RANDOM_NT, placeholder=True)]


SEQ_PRIMERS = (
    sp.custom("Read 1", "TruSeq Read 1", il.TRUSEQ_READ1,
              "Illumina adapter sequence"),
    sp.custom("Index 1 (i7)", "TruSeq Index 1", il.INDEX1_PRIMER,
              "Illumina adapter sequence"),
    sp.custom("Read 2", "TruSeq Read 2", il.TRUSEQ_READ2,
              "Illumina adapter sequence"),
)


def nome_library() -> Construct:
    """Outer published primer arms with one explicit unknown kit-derived junction each."""
    con = Construct([
        seg("P5", il.P5, "p5"), seg("Read 1 remainder", il.TRUSEQ_READ1[4:], "r1"),
        seg("kit remnant", "?", inferred=True, placeholder=True,
            note="Pico Methyl-Seq internal primer sequence and length not published"),
        seg("bisulfite insert", "XXXXXXXX...XXXXXXXX", placeholder=True),
        seg("kit remnant'", "?", inferred=True, placeholder=True,
            note="Pico Methyl-Seq internal primer sequence and length not published"),
        seg("Read 2 arm'", revcomp(il.TRUSEQ_READ2), "r2"),
        seg("i7'", "I" * INDEX_NT_NOME, "cbc", placeholder=True,
            feature=feature("sample_i7", "sample_index", "unknown")),
        seg("P7'", il.P7_RC, "p7"),
    ], name="scNOMe-seq library")
    _validate(con)
    return con


def cool_library() -> Construct:
    arm = revcomp(il.TRUSEQ_READ2)
    known = revcomp(COOL_P2_HANDLE)
    if not arm.startswith(known):
        raise ValueError("primer-2 handle is not a suffix of TruSeq Read 2")
    con = Construct([
        seg("P5", il.P5, "p5", inferred=True),
        seg("Read 1 remainder", il.TRUSEQ_READ1[4:], "r1", inferred=True),
        seg("primer-1 N9", "N" * RANDOM_NT, placeholder=True),
        seg("bisulfite insert", "XXXXXXXX...XXXXXXXX", placeholder=True),
        seg("primer-2 N9'", "N" * RANDOM_NT, placeholder=True),
        seg("primer-2 handle'", known, "r2"),
        seg("reverse-primer extension", arm[len(known):], "r2", inferred=True,
            note="standard TruSeq-shaped NEB indexed primer; exact product not printed"),
        seg("i7'", "I" * INDEX_NT_COOL, "cbc", placeholder=True, inferred=True,
            feature=feature("sample_i7", "sample_index", "unknown")),
        seg("P7'", il.P7_RC, "p7", inferred=True),
    ], name="scCOOL-seq library")
    _validate(con)
    return con


def _validate(con: Construct) -> None:
    for primer in SEQ_PRIMERS:
        if sp.locate(con, primer) is None:
            raise ValueError(f"{con.name}: {primer.name} site is absent")


def library(protocol: str) -> Construct:
    if protocol == SCNOME: return nome_library()
    if protocol == SCCOOL: return cool_library()
    raise ValueError(protocol)
