"""Molecular constructors for SCRB-seq and mcSCRB-seq."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Segment, feature
import illumina as il
import nextera as nx
import seqprimers as sp

SCRB = "SCRB-seq"
MCSCRB = "mcSCRB-seq"
PROTOCOLS = (SCRB, MCSCRB)

FULL_HANDLE = il.TRUSEQ_READ1
PCR_HANDLE = FULL_HANDLE[:22]
BARCODE_LEN = 6
UMI_LEN = 10
DT_LEN = 30
P5_ENRICH = il.P5 + FULL_HANDLE[4:]
BLOCKED_TSO_BASES = "CGC" + PCR_HANDLE + "GGG"
UNBLOCKED_TSO_BASES = PCR_HANDLE + "GGG"

SOURCE_E3 = "ACACTCTTTCCCTACACGACGCTCTTCCGATCTBBBBBBNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN"
SOURCE_MC_TSO = "ACACTCTTTCCCTACACGACGCGGG"
CELL_BARCODE = feature("cell_barcode", "cell_barcode", "whitelist",
                       whitelist="scrb-seq-e3-primer-set")
UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "unknown")
FEATURES = {"cell barcode": CELL_BARCODE, "UMI": UMI_FEATURE, "i7'": I7_FEATURE}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


def oligo_dt_segments() -> list[Segment]:
    return [seg("TruSeq Read 1", FULL_HANDLE, "r1"),
            seg("cell barcode", "B" * BARCODE_LEN, "cbc", placeholder=True),
            seg("UMI", "N" * UMI_LEN, "umi", placeholder=True),
            seg("T30", "T" * DT_LEN), seg("VN", "VN", placeholder=True)]


def tso_segments(protocol: str) -> list[Segment]:
    if protocol == SCRB:
        return [seg("iso-base block", "iCiGiC", placeholder=True),
                seg("PCR handle", PCR_HANDLE, "r1"),
                seg("rGrGrG", "rGrGrG", placeholder=True)]
    if protocol == MCSCRB:
        return [seg("PCR handle", PCR_HANDLE, "r1"),
                seg("rGrGrG", "rGrGrG", placeholder=True)]
    raise ValueError(protocol)


READ_PRIMERS = (
    sp.custom("Read 1", "TruSeq Read 1", il.TRUSEQ_READ1,
              "Illumina adapter sequence; site installed by E3V6NEXT/P5NEXTPT5"),
    sp.NEXTERA["I1"],
    sp.NEXTERA["R2"],
)


def final_library(protocol: str) -> Construct:
    if protocol not in PROTOCOLS:
        raise ValueError(protocol)
    con = Construct([
        seg("P5", il.P5, "p5"),
        seg("Read 1 remainder", FULL_HANDLE[4:], "r1"),
        seg("cell barcode", "B" * BARCODE_LEN, "cbc", placeholder=True),
        seg("UMI", "N" * UMI_LEN, "umi", placeholder=True),
        seg("poly(T)", "T" * DT_LEN), seg("VN", "VN", placeholder=True),
        seg("cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
        seg("ME'", nx.ME_RC, "me"), seg("s7'", nx.S7_RC, "s7"),
        seg("i7'", "I" * 8, "cbc", placeholder=True), seg("P7'", il.P7_RC, "p7"),
    ], name=f"{protocol} library")
    for primer in READ_PRIMERS:
        if sp.locate(con, primer) is None:
            raise ValueError(f"{protocol}: {primer.name} site is absent")
    return con


def assembled_e3() -> str:
    return "".join(s.top for s in oligo_dt_segments())

