"""Construct model shared by the individual CEL-Seq and CEL-Seq2 pages.

The primer designs are transcribed from Table S2 of the CEL-Seq2 paper.  That table is
also the accessible primary source for the original CEL-Seq primer design; exact TruSeq
Small RNA library oligos are not present in either accessible primary-source package and
are deliberately not supplied here.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, complement_segments, feature

CELL_BARCODE = feature("cell_barcode", "cell_barcode", "whitelist",
                       whitelist="cel-seq-rt-primer-set")
UMI_FEATURE = feature("umi", "umi", "random")
FEATURES = {"cell barcode": CELL_BARCODE, "UMI": UMI_FEATURE}

T7_PROMOTER = "TAATACGACTCACTATAGGG"
T7_TRANSCRIBED_G = "GGG"
POLYT_NT = 24
ANCHOR = "V"

CEL1_SPACER = "CGATTGAGGCCGG"
CEL1_RA5 = "GTTCAGAGTTCTACAGTCCGACGATC"
CEL1_BARCODE_NT = 8

CEL2_SPACER = "GCCGG"
CEL2_RA5 = "AGTTCTACAGTCCGACGATC"
CEL2_UMI_NT = 6
CEL2_BARCODE_NT = 6

RANDOMHEX_TAIL = "GCCTTGGCACCCGAGAATTCCA"
RANDOMHEX_NT = 6
RANDOMHEX_RT = RANDOMHEX_TAIL + "N" * RANDOMHEX_NT


def _seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


def rt_primer_segments(protocol: str) -> list[Segment]:
    """Published primer design with invalid variants made unrepresentable."""
    if protocol == "CEL-Seq":
        spacer, ra5, umi_n, barcode_n = CEL1_SPACER, CEL1_RA5, 0, CEL1_BARCODE_NT
    elif protocol == "CEL-Seq2":
        spacer, ra5, umi_n, barcode_n = CEL2_SPACER, CEL2_RA5, CEL2_UMI_NT, CEL2_BARCODE_NT
    else:
        raise ValueError(f"unknown CEL-Seq protocol: {protocol!r}")
    out = [
        _seg("5' spacer", spacer),
        _seg("T7 promoter", T7_PROMOTER, "t7"),
        _seg("small-RNA 5' sequence", ra5, "r1"),
    ]
    if umi_n:
        out.append(_seg("UMI", "N" * umi_n, "umi", placeholder=True))
    out.extend([
        _seg("cell barcode", "N" * barcode_n, "cbc", placeholder=True),
        _seg("oligo-dT", "T" * POLYT_NT),
        _seg("V anchor", ANCHOR, placeholder=True),
    ])
    return out


def rt_primer(protocol: str) -> Construct:
    return Construct(rt_primer_segments(protocol), name=f"{protocol} RT primer")


def rt_scene(protocol: str) -> Scene:
    mrna = [_seg("mRNA body", "XXXXXXXXXXXX", placeholder=True),
            _seg("poly(A) junction", "B", placeholder=True),
            _seg("poly(A)", "A" * POLYT_NT)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA")
    sc.anneal("RT primer", rt_primer_segments(protocol), to="mRNA",
              pair=("oligo-dT", "poly(A)"), label="RT primer")
    sc.mark("RT primer", "V anchor", "anchored at the poly(A) junction")
    sc.arrow("RT primer", "reverse transcriptase")
    return sc


def double_stranded_cdna(protocol: str) -> Construct:
    return Construct([
        *rt_primer_segments(protocol),
        _seg("cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
    ], name=f"{protocol} double-stranded cDNA")


def arna(protocol: str) -> Construct:
    """The amplified antisense RNA, downstream of the double-stranded T7 promoter."""
    primer = rt_primer_segments(protocol)
    first_ra5 = next(i for i, s in enumerate(primer) if s.name == "small-RNA 5' sequence")
    transcribed = primer[first_ra5:]
    return Construct([
        _seg("T7 initial G", T7_TRANSCRIBED_G, "t7"),
        *transcribed,
        _seg("antisense RNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
    ], name=f"{protocol} amplified RNA")


def cel2_random_rt_scene() -> Scene:
    """Random-hexamer priming on the fragmented, 5'-most CEL-Seq2 aRNA piece."""
    fragment = [
        _seg("5' tag", "GGG" + CEL2_RA5, "r1"),
        _seg("UMI", "N" * CEL2_UMI_NT, "umi", placeholder=True),
        _seg("cell barcode", "N" * CEL2_BARCODE_NT, "cbc", placeholder=True),
        _seg("poly(U)/insert", "UUUXXXXXXXXXXXX", placeholder=True),
        _seg("random site", "XXXXXX", placeholder=True),
    ]
    primer = [_seg("library-RT tail", RANDOMHEX_TAIL, "r2"),
              _seg("random hexamer", "NNNNNN", placeholder=True)]
    sc = Scene()
    sc.strand("aRNA fragment", fragment, label="aRNA")
    sc.anneal("randomhexRT", primer, to="aRNA fragment",
              pair=("random hexamer", "random site"), label="randomhexRT")
    sc.arrow("randomhexRT", "reverse transcriptase")
    return sc


def cel1_ligated_arna() -> Construct:
    """Authoritative boundary: original aRNA ligation is known, adaptor bases are not."""
    return Construct([
        _seg("5' tag", T7_TRANSCRIBED_G + CEL1_RA5, "r1"),
        _seg("cell barcode", "N" * CEL1_BARCODE_NT, "cbc", placeholder=True),
        _seg("poly(U)/insert", "UUUXXXXXXXX...XXXXXXXX", placeholder=True),
        _seg("unpublished RNA 3' adaptor", "[RNA 3' ADAPTOR]", placeholder=True,
             inferred=True),
    ], name="CEL-Seq ligated aRNA")


def final_library(protocol: str) -> Construct:
    """Supported inner architecture plus explicitly unpublished outer PCR arms."""
    if protocol == "CEL-Seq":
        middle = [
            _seg("small-RNA 5' sequence", CEL1_RA5, "r1"),
            _seg("cell barcode", "N" * CEL1_BARCODE_NT, "cbc", placeholder=True),
            _seg("poly(T)/cDNA", "TTTXXXXXXXX...XXXXXXXX", placeholder=True),
            _seg("unpublished ligated-adaptor sequence", "[3' ADAPTOR]", inferred=True,
                 placeholder=True),
        ]
    elif protocol == "CEL-Seq2":
        middle = [
            _seg("RT-primer 5' sequence", CEL2_RA5, "r1"),
            _seg("UMI", "N" * CEL2_UMI_NT, "umi", placeholder=True),
            _seg("cell barcode", "N" * CEL2_BARCODE_NT, "cbc", placeholder=True),
            _seg("poly(T)/cDNA", "TTTXXXXXXXX...XXXXXXXX", placeholder=True),
            _seg("randomhexRT tail", RANDOMHEX_TAIL, "r2"),
        ]
    else:
        raise ValueError(f"unknown CEL-Seq protocol: {protocol!r}")
    return Construct([
        _seg("unpublished P5-side PCR arm", "[P5-SIDE ARM]", inferred=True,
             placeholder=True),
        *middle,
        _seg("unpublished indexed P7-side PCR arm", "[INDEXED P7-SIDE ARM]",
             inferred=True, placeholder=True),
    ], name=f"{protocol} final library, authoritative boundary")


def _validate_model() -> None:
    cel1 = rt_primer("CEL-Seq")
    cel2 = rt_primer("CEL-Seq2")
    expected1 = (CEL1_SPACER + T7_PROMOTER + CEL1_RA5 + "N" * CEL1_BARCODE_NT
                 + "T" * POLYT_NT + ANCHOR)
    expected2 = (CEL2_SPACER + T7_PROMOTER + CEL2_RA5 + "N" * CEL2_UMI_NT
                 + "N" * CEL2_BARCODE_NT + "T" * POLYT_NT + ANCHOR)
    if cel1.top() != expected1 or cel2.top() != expected2:
        raise ValueError("RT primer assembly no longer matches the source designs")
    if RANDOMHEX_RT != RANDOMHEX_TAIL + "N" * RANDOMHEX_NT:
        raise ValueError("randomhexRT tail/hexamer boundary is invalid")
    rt_scene("CEL-Seq").rows()
    rt_scene("CEL-Seq2").rows()
    cel2_random_rt_scene().rows()


_validate_model()
