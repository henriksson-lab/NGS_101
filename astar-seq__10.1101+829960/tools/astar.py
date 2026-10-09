"""Sequence model for ASTAR-seq.

The C1 oligos and custom ATAC indexing-primer structures are from Supplementary
Table 5 of Xing et al. (bioRxiv 829960).  The RNA library uses the same Nextera
back end as the ATAC library after the separated cDNA is tagmented.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Segment, feature, revcomp
from illumina import P5, P7
import nextera as nx
import seqprimers as sp


C1_HANDLE = "GGCGACAACACCGATTGATCA"
C1_DT_SPACER = "CG"
C1_DT_LEN = 31
C1_DT = C1_HANDLE + C1_DT_SPACER + "T" * C1_DT_LEN
C1_TSO = C1_HANDLE + "GGG"       # ordered as an all-RNA oligo
C1_PCR = C1_HANDLE

# Modifications as printed in Supplementary Table 5.
C1_DT_MOD = "/5BiosG/"
C1_PCR_MOD = "/5BiosG/"
C1_TSO_MOD = "all bases are ribonucleotides"

ATAC_QPCR_S5 = nx.ADAPTOR_S5
ATAC_QPCR_S7 = nx.ADAPTOR_S7
INDEX_LEN = 8


def _seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def c1_cdna() -> Construct:
    """Representative full-length cDNA after template switching and C1 PCR.

    The transcript body is a placeholder; both terminal handles derive from the
    shared sequence on the oligo-dT and TSO.
    """
    return Construct([
        _seg("C1 handle, oligo-dT end", C1_HANDLE, "tso"),
        _seg("CG", C1_DT_SPACER),
        _seg("poly(T)", "T" * C1_DT_LEN),
        _seg("cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
        _seg("untemplated CCC", "CCC", "tso"),
        _seg("C1 handle, TSO end", revcomp(C1_HANDLE), "tso"),
    ], name="ASTAR full-length cDNA")


def nextera_library(insert_name: str = "insert") -> Construct:
    """Dual-indexed Nextera library used for either ASTAR modality."""
    lib = Construct([
        _seg("P5", P5, "p5"),
        _seg("i5", "N" * INDEX_LEN, "cbc", placeholder=True,
             feature=feature("sample_i5", "sample_index", "unknown")),
        _seg("s5", nx.S5, "s5"),
        _seg("ME", nx.ME, "me"),
        _seg(insert_name, "XXXXXXXX...XXXXXXXX", placeholder=True),
        _seg("ME'", nx.ME_RC, "me"),
        _seg("s7'", nx.S7_RC, "s7"),
        _seg("i7'", "N" * INDEX_LEN, "cbc", placeholder=True,
             feature=feature("sample_i7", "sample_index", "unknown")),
        _seg("P7'", revcomp(P7), "p7"),
    ], name=f"ASTAR {insert_name} library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError(f"invalid ASTAR {insert_name} library: " + "; ".join(problems))
    return lib


def atac_library() -> Construct:
    return nextera_library("accessible genomic DNA")


def rna_library() -> Construct:
    return nextera_library("cDNA")


SEQ_PRIMERS = tuple(sp.NEXTERA[k] for k in ("R1", "I1", "I2", "R2"))
