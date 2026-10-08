"""Sequence model for Smallwood et al. scBS-seq (Nature Methods, 2014).

The defining paper prints oligo1, oligo2 and PE1.0 and cites Quail et al. for the
indexed iPCRTag primer.  Quail's Supplementary Table 1 prints that primer family,
the PE adapter and the dedicated index-read primer.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Segment, revcomp


def _seg(name, top, tag=None, **kw):
    return Segment(name=name, top=top, tag=tag, **kw)


# Printed verbatim in the scBS-seq Online Methods.
OLIGO1_HANDLE = "CTACACGACGCTCTTCCGATCT"
OLIGO2_HANDLE = "TGCTGAACCGCTCTTCCGATCT"
RANDOM_NT = 9
PE1 = il.P5 + il.TRUSEQ_READ1[4:]  # the shared ACAC is present only once

# Quail et al. 2012, doi:10.1038/nmeth.1814, Supplementary Table 1.
IPCRTAG_INDEX_NT = 8
READ2_PRIMER = "CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"
IPCRTAG_3_ARM = "GAGAT" + READ2_PRIMER
INDEX_PRIMER = "AAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTC"


def oligo1_segments() -> list[Segment]:
    return [_seg("oligo1 handle", OLIGO1_HANDLE, "r1"),
            _seg("N9", "N" * RANDOM_NT, placeholder=True)]


def oligo2_segments() -> list[Segment]:
    return [_seg("oligo2 handle", OLIGO2_HANDLE, "r2"),
            _seg("N9", "N" * RANDOM_NT, placeholder=True)]


def ipcrtag_primer_segments() -> list[Segment]:
    """Quail iPCRTag primer family; the eight-base member index is variable."""
    return [
        _seg("P7", il.P7, "p7"),
        _seg("i7", "I" * IPCRTAG_INDEX_NT, "cbc", placeholder=True),
        _seg("3' arm", IPCRTAG_3_ARM, "r2"),
    ]


def final_library(insert_nt: int = 28) -> Construct:
    """Indexed library, top strand 5'->3' from P5 to P7'.

    A short placeholder insert keeps the page readable.  The two N9s are random
    priming sequence, not UMIs, despite sharing the UMI colour token.
    """
    insert = "X" * insert_nt if insert_nt <= 20 else "XXXXXXXXXX...XXXXXXXXXX"
    lib = Construct([
        # P5 and Read 1 overlap by ACAC in PE1, so segment the 58-nt primer this way.
        _seg("P5", il.P5, "p5"),
        _seg("Read 1 remainder", il.TRUSEQ_READ1[4:], "r1"),
        _seg("oligo1 N9", "N" * RANDOM_NT, placeholder=True),
        _seg("bisulfite insert", insert, placeholder=True),
        _seg("oligo2 N9", "N" * RANDOM_NT, placeholder=True),
        _seg("index-read site", revcomp(IPCRTAG_3_ARM), "r2"),
        _seg("i7", "I" * IPCRTAG_INDEX_NT, "cbc", placeholder=True),
        _seg("P7'", il.P7_RC, "p7"),
    ], name="scBS-seq indexed library")
    problems = sp.verify(lib, SEQ_PRIMERS,
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if problems:
        raise ValueError("invalid scBS-seq final library: " + "; ".join(problems))
    return lib


SEQ_PRIMERS = [
    sp.custom("Read 1", "TruSeq Read 1", il.TRUSEQ_READ1,
              'Illumina "Illumina Adapter Sequences" #1000000002694'),
    sp.custom("Index 1 (i7)", "historical iPCRTag index primer", INDEX_PRIMER,
              "Quail et al. 2012 Supplementary Table 1"),
    sp.custom("Read 2", "historical paired-end Read 2", READ2_PRIMER,
              'Illumina "Oligonucleotide Sequences for Paired End DNA" (obsolete)'),
]
