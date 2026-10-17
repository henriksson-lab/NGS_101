"""Construct model shared by the individual CEL-Seq and CEL-Seq2 pages.

The CEL-Seq2 designs are transcribed from Table S2 of the CEL-Seq2 paper.  The original
CEL-Seq design (13-nt spacer, T7 promoter, RA5, 8-nt barcode, T24V) is printed for all 40
primers in the 2012 detailed "CEL-Seq Protocol" (Cell Reports Supplemental Information,
mmc4) and agrees with Table S2.  That protocol names the Illumina TruSeq Small RNA kit
reagents (RA3, RTP, RP1, RPIX) but does not print their bases, so exact library oligos are
deliberately not supplied here.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from adapter_ligation import AdapterLigationGeometry
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

@dataclass(frozen=True)
class PrimerDesign:
    spacer: str
    ra5: str
    umi_nt: int
    barcode_nt: int
    randomhex_tail: str = ""
    randomhex_nt: int = 0


CEL1_DESIGN = PrimerDesign(CEL1_SPACER, CEL1_RA5, 0, CEL1_BARCODE_NT)


def _seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


def rt_primer_segments(protocol: str, design: PrimerDesign | None = None) -> list[Segment]:
    """Published primer design with invalid variants made unrepresentable."""
    design = design or (CEL1_DESIGN if protocol == "CEL-Seq" else None)
    if design is None:
        raise ValueError(f"{protocol}: supply the paper's primer design")
    spacer, ra5, umi_n, barcode_n = (design.spacer, design.ra5,
                                     design.umi_nt, design.barcode_nt)
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


def rt_primer(protocol: str, design: PrimerDesign | None = None) -> Construct:
    return Construct(rt_primer_segments(protocol, design), name=f"{protocol} RT primer")


def rt_scene(protocol: str, design: PrimerDesign | None = None) -> Scene:
    mrna = [_seg("mRNA body", "XXXXXXXXXXXX", placeholder=True),
            _seg("poly(A) junction", "B", placeholder=True),
            _seg("poly(A)", "A" * POLYT_NT)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA")
    sc.anneal("RT primer", rt_primer_segments(protocol, design), to="mRNA",
              pair=("oligo-dT", "poly(A)"), label="RT primer")
    sc.mark("RT primer", "V anchor", "anchored at the poly(A) junction")
    sc.arrow("RT primer", "reverse transcriptase")
    return sc


def double_stranded_cdna(protocol: str, design: PrimerDesign | None = None) -> Construct:
    return Construct([
        *rt_primer_segments(protocol, design),
        _seg("cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
    ], name=f"{protocol} double-stranded cDNA")


def arna(protocol: str, design: PrimerDesign | None = None) -> Construct:
    """The amplified antisense RNA, downstream of the double-stranded T7 promoter."""
    primer = rt_primer_segments(protocol, design)
    first_ra5 = next(i for i, s in enumerate(primer) if s.name == "small-RNA 5' sequence")
    transcribed = primer[first_ra5:]
    return Construct([
        _seg("T7 initial G", T7_TRANSCRIBED_G, "t7"),
        *transcribed,
        _seg("antisense RNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
    ], name=f"{protocol} amplified RNA")


def random_rt_scene(design: PrimerDesign) -> Scene:
    """Random-hexamer priming on the fragmented, 5'-most CEL-Seq2 aRNA piece."""
    fragment = [
        _seg("5' tag", "GGG" + design.ra5, "r1"),
        _seg("UMI", "N" * design.umi_nt, "umi", placeholder=True),
        _seg("cell barcode", "N" * design.barcode_nt, "cbc", placeholder=True),
        _seg("poly(U)/insert", "UUUXXXXXXXXXXXX", placeholder=True),
        _seg("random site", "X" * design.randomhex_nt, placeholder=True),
    ]
    primer = [_seg("library-RT tail", design.randomhex_tail, "r2"),
              _seg("random hexamer", "N" * design.randomhex_nt, placeholder=True)]
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


def cel1_adapter_geometry() -> AdapterLigationGeometry:
    return AdapterLigationGeometry.unavailable(
        "The 2012 CEL-Seq detailed protocol ligates the Illumina TruSeq Small RNA kit 3′ "
        "adapter (RA3, diluted 5-fold) to phosphatase- and PNK-treated fragmented aRNA with "
        "T4 RNA ligase 2, truncated (28 °C, 1 h), but prints neither the RA3 bases nor its "
        "terminal modifications or ligation-facing end geometry (\"sequences available from "
        "Illumina\").")


def final_library(protocol: str, design: PrimerDesign | None = None) -> Construct:
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
        if design is None:
            raise ValueError("CEL-Seq2 requires its source-specific primer design")
        middle = [
            _seg("RT-primer 5' sequence", design.ra5, "r1"),
            _seg("UMI", "N" * design.umi_nt, "umi", placeholder=True),
            _seg("cell barcode", "N" * design.barcode_nt, "cbc", placeholder=True),
            _seg("poly(T)/cDNA", "TTTXXXXXXXX...XXXXXXXX", placeholder=True),
            _seg("randomhexRT tail", design.randomhex_tail, "r2"),
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
    expected1 = (CEL1_SPACER + T7_PROMOTER + CEL1_RA5 + "N" * CEL1_BARCODE_NT
                 + "T" * POLYT_NT + ANCHOR)
    if cel1.top() != expected1:
        raise ValueError("RT primer assembly no longer matches the source designs")
    rt_scene("CEL-Seq").rows()


_validate_model()
