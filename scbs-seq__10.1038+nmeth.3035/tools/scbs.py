"""Sequence model for Smallwood et al. scBS-seq (Nature Methods, 2014).

The defining paper prints oligo1, oligo2 and PE1.0.  It names the indexed iPCRTag
primer but refers to Quail et al. for its sequence; that historical arm is therefore
marked inferred wherever it is used in the final construct.
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

# Secondary sequence source: the upstream scg_lib_structs rendering of the Quail
# iPCRTag design.  Keep it separate from the paper-derived constants above.
IPCRTAG_INDEX_NT = 8
READ2_PRIMER = "CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"
IPCRTAG_3_ARM = "GAGAT" + READ2_PRIMER
INDEX_PRIMER = revcomp(IPCRTAG_3_ARM)


def oligo1_segments() -> list[Segment]:
    return [_seg("oligo1 handle", OLIGO1_HANDLE, "r1"),
            _seg("N9", "N" * RANDOM_NT, placeholder=True)]


def oligo2_segments() -> list[Segment]:
    return [_seg("oligo2 handle", OLIGO2_HANDLE, "r2"),
            _seg("N9", "N" * RANDOM_NT, placeholder=True)]


def ipcrtag_primer_segments() -> list[Segment]:
    """Secondary-source primer; inference styling is inseparable from construction."""
    return [
        _seg("P7", il.P7, "p7", inferred=True),
        _seg("i7", "I" * IPCRTAG_INDEX_NT, "cbc", placeholder=True, inferred=True),
        _seg("3' arm", IPCRTAG_3_ARM, "r2", inferred=True),
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
        _seg("index-read site", INDEX_PRIMER, "r2", inferred=True,
             note="iPCRTag arm from the secondary upstream page"),
        _seg("i7", "I" * IPCRTAG_INDEX_NT, "cbc", placeholder=True, inferred=True,
             note="indexed iPCRTag named by the paper; index sequence not printed"),
        _seg("P7'", il.P7_RC, "p7", inferred=True,
             note="iPCRTag sequence supplied by the secondary upstream page"),
    ], name="scBS-seq indexed library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid scBS-seq final library: " + "; ".join(problems))
    return lib


SEQ_PRIMERS = [
    sp.custom("Read 1", "TruSeq Read 1", il.TRUSEQ_READ1,
              'Illumina "Illumina Adapter Sequences" #1000000002694'),
    sp.custom("Index 1 (i7)", "historical iPCRTag index primer", INDEX_PRIMER,
              "Quail et al. design, via upstream scg_lib_structs"),
    sp.SeqPrimer("Index 2 (i5)", "TruSeq Index 2 (no i5 index)",
                 il.INDEX2_PRIMER_RC,
                 'Illumina "Indexed Sequencing Overview Guide" #15057455',
                 "The primer can bind, but this single-indexed library has no i5 index."),
    sp.custom("Read 2", "historical paired-end Read 2", READ2_PRIMER,
              "Quail et al. design, via upstream scg_lib_structs"),
]
