"""Molecular model for snATAC-seq (Preissl et al. 2018)."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Segment, feature, revcomp
import illumina as il
import nextera as nx
import seqprimers as sp


# Constant portions of the two barcoded transposons, Supplementary Table 5.
P5_LINKER = "TCCACGC"
P7_LINKER = "CTGTCCCTGTCC"
READ1_SITE = "GCGATCGAGGACGGC"
READ2_SITE = "CACCGTCTCCGCCTC"
BARCODE_LEN = 8
P5_BARCODES = 8
P7_BARCODES = 12
PCR_INDEX_LEN = 8

PMENTS = nx.ME_RC
READ1_PRIMER = READ1_SITE + nx.ME
READ2_PRIMER = READ2_SITE + nx.ME
INDEX1_PRIMER = revcomp(READ2_PRIMER)

# Representative verbatim rows retained for the small transcription test.
P5_1 = "TCGTCGGCAGCGTCTCCACGCTATAGCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG"
P7_1 = "GTCTCGTGGGCTCGGCTGTCCCTGTCCCGAGTAATCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG"
N701 = "CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGG"
S502 = "AATGATACGGCGACCACCGAGATCTACACCTCTCTATTCGTCGGCAGCGTC"
P5_1_BARCODE = "TATAGCCT"
P7_1_BARCODE = "CGAGTAAT"
N701_ORDERED_INDEX = "TCGCCTTA"
S502_INDEX = "CTCTCTAT"


def _seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def barcode(name: str, value: str | None, letter: str) -> Segment:
    """A real or placeholder barcode whose required length is enforced here."""
    seq = letter * BARCODE_LEN if value is None else value
    if len(seq) != BARCODE_LEN:
        raise ValueError(f"{name} must be {BARCODE_LEN} nt")
    role = "cell_barcode" if name in ("p5 barcode", "p7 barcode") else "cell_barcode"
    part = name.replace(" barcode", "")
    return _seg(name, seq, "cbc", placeholder=value is None,
                feature=feature(f"cell_{part.replace(' ', '_')}", role, "combinatorial",
                                group="cell_id", part=part))


def p5_transposon(value: str | None = None) -> list[Segment]:
    return [_seg("s5", nx.S5, "s5"), _seg("p5 linker", P5_LINKER, "r2"),
            barcode("p5 barcode", value, "A"), _seg("Read 1 site", READ1_SITE, "r1"),
            _seg("ME", nx.ME, "me")]


def p7_transposon(value: str | None = None) -> list[Segment]:
    return [_seg("s7", nx.S7, "s7"), _seg("p7 linker", P7_LINKER, "r3"),
            barcode("p7 barcode", value, "B"), _seg("Read 2 site", READ2_SITE, "r1"),
            _seg("ME", nx.ME, "me")]


def tagged_fragment() -> Construct:
    """The amplifiable A...B tagmentation product after its 9-nt gaps are filled."""
    right = p7_transposon()
    return Construct([
        *p5_transposon(),
        _seg("genomic insert", "XXXXXXXX...XXXXXXXX", placeholder=True),
        *[Segment(s.name + "'", revcomp(s.top) if not s.placeholder else s.top.lower(),
                  s.tag, s.placeholder, feature=s.feature) for s in reversed(right)],
    ], name="gap-filled snATAC fragment")


SEQ_PRIMERS = (
    sp.custom("Read 1", "snATAC Read 1", READ1_PRIMER,
              "Preissl et al. Supplementary Table 5"),
    sp.custom("Index 1 (i7)", "snATAC Index 1", INDEX1_PRIMER,
              "Preissl et al. Supplementary Table 5"),
    sp.custom("Index 2 (i5)", "flow-cell P5 oligo", il.P5,
              "HiSeq 2500 forward-strand index workflow",
              "No separate Index 2 sequencing primer is listed."),
    sp.custom("Read 2", "snATAC Read 2", READ2_PRIMER,
              "Preissl et al. Supplementary Table 5"),
)


def final_library() -> Construct:
    """Dual-PCR-indexed, dual-transposon-indexed snATAC library.

    The B-side sequence is generated from the ordered p7 oligo rather than retyped.
    Primer placement is validated here so no caller can obtain an unsequenceable model.
    """
    tagged = tagged_fragment()
    con = Construct([
        _seg("P5", il.P5, "p5"),
        _seg("i5", "N" * PCR_INDEX_LEN, "cbc", placeholder=True,
             feature=feature("cell_i5", "cell_barcode", "combinatorial", group="cell_id", part="i5")),
        *tagged.segments,
        _seg("i7 read", "N" * PCR_INDEX_LEN, "cbc", placeholder=True,
             feature=feature("cell_i7", "cell_barcode", "combinatorial", group="cell_id", part="i7")),
        _seg("P7'", il.P7_RC, "p7"),
    ], name="snATAC-seq library")
    errors = sp.verify(con, SEQ_PRIMERS)
    if errors:
        raise ValueError("invalid snATAC library: " + "; ".join(errors))
    return con


def p5_index_primer(value: str | None = None) -> list[Segment]:
    idx = barcode("i5", value, "N")
    return [_seg("P5", il.P5, "p5"), idx,
            _seg("s5", nx.S5, "s5")]


def p7_index_primer(value: str | None = None) -> list[Segment]:
    idx = barcode("i7 as ordered", value, "N")
    return [_seg("P7", il.P7, "p7"),
            idx,
            _seg("s7", nx.S7, "s7")]
