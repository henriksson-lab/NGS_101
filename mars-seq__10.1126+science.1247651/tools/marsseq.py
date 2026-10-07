"""Construct model for the individual MARS-seq protocol pages."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, revcomp

T7 = "TAATACGACTCACTATAGGG"
SPACER = "CGATTGAGGCCGG"
MARS2_TRANSCRIBED_LEADER = "GC"
MARS2_RT_HANDLE = "GACGTGTGCTCTTCCGATCT"
MARS2_WELL_NT = 7
MARS2_UMI_NT = 8
POLYT_NT = 20

POOL_ORDERED = "GACT"                 # representative ix1, as ordered
POOL_READ = revcomp(POOL_ORDERED)      # AGTC, Supplementary Table 2 "Seq barcode"
POOL_RANDOM_NT = 5
LIG_CONSTANT = "AGATCGGAAGAGCGTCGTGTAG"


def _seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def rt1_segments(protocol: str) -> list[Segment]:
    if protocol == "MARS-seq2.0":
        return [
            _seg("5' spacer", SPACER), _seg("T7 promoter", T7, "t7"),
            _seg("transcribed leader", MARS2_TRANSCRIBED_LEADER),
            _seg("RT/PCR handle", MARS2_RT_HANDLE, "r2"),
            _seg("well barcode", "N" * MARS2_WELL_NT, "cbc", placeholder=True),
            _seg("UMI", "N" * MARS2_UMI_NT, "umi", placeholder=True),
            _seg("oligo-dT", "T" * POLYT_NT),
            _seg("N anchor", "N", placeholder=True),
        ]
    if protocol == "MARS-seq":
        # The 2014 supplement containing the bases was unavailable.  Keep only the
        # architecture established by the paper, and do not copy the upstream sequence.
        return [
            _seg("unavailable T7/RT handle", "XXXXXXXXXXXXXX", inferred=True,
                 placeholder=True),
            _seg("well barcode", "XXXXXX", "cbc", inferred=True, placeholder=True),
            _seg("UMI", "XXXX", "umi", inferred=True, placeholder=True),
            _seg("oligo-dT", "T" * POLYT_NT, inferred=True),
            _seg("anchor", "N", inferred=True, placeholder=True),
        ]
    raise ValueError(f"unknown MARS-seq protocol: {protocol!r}")


def rt_scene(protocol: str) -> Scene:
    mrna = [_seg("body", "XXXXXXXXXXXX", placeholder=True),
            _seg("junction", "B", placeholder=True),
            _seg("poly(A)", "A" * POLYT_NT)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA")
    sc.anneal("RT1", rt1_segments(protocol), to="mRNA",
              pair=("oligo-dT", "poly(A)"), label="RT1")
    sc.mark("RT1", "anchor" if protocol == "MARS-seq" else "N anchor",
            "poly(A) junction")
    sc.arrow("RT1", "reverse transcriptase")
    return sc


def ds_cdna(protocol: str) -> Construct:
    return Construct([*rt1_segments(protocol),
                      _seg("cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True)],
                     name=f"{protocol} double-stranded cDNA")


def mars2_ligation_adapter() -> Construct:
    return Construct([
        _seg("ordered pool barcode", POOL_ORDERED, "cbc"),
        _seg("random diversity", "N" * POOL_RANDOM_NT, placeholder=True),
        _seg("ligation/RT2 handle", LIG_CONSTANT, "r1"),
    ], name="MARS-seq2.0 ligation adapter")


def ligated_arna(protocol: str) -> Construct:
    if protocol == "MARS-seq2.0":
        adaptor = list(mars2_ligation_adapter())
    elif protocol == "MARS-seq":
        adaptor = [_seg("unavailable pool-barcoded adapter", "[POOL ADAPTOR]",
                        inferred=True, placeholder=True)]
    else:
        raise ValueError(f"unknown MARS-seq protocol: {protocol!r}")
    return Construct([
        _seg("barcoded aRNA 5' end", "[WELL BC / UMI / poly(U)]", placeholder=True,
             inferred=protocol == "MARS-seq"),
        _seg("antisense insert", "XXXXXXXX...XXXXXXXX", placeholder=True),
        *adaptor,
    ], name=f"{protocol} ligated aRNA")


def final_library(protocol: str) -> Construct:
    if protocol == "MARS-seq2.0":
        middle = [
            _seg("random diversity", "N" * POOL_RANDOM_NT, placeholder=True),
            _seg("pool barcode", POOL_READ, "cbc"),
            _seg("sense cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
            _seg("poly(A)", "AAA", placeholder=True),
            _seg("UMI'", "N" * MARS2_UMI_NT, "umi", placeholder=True),
            _seg("well barcode'", "N" * MARS2_WELL_NT, "cbc", placeholder=True),
        ]
    elif protocol == "MARS-seq":
        middle = [_seg("pool barcode / cDNA / UMI / well barcode",
                       "[POOL BC]XXXXXXXX...[UMI][CELL BC]", inferred=True,
                       placeholder=True)]
    else:
        raise ValueError(f"unknown MARS-seq protocol: {protocol!r}")
    return Construct([
        _seg("unpublished P5-side PCR arm", "[P5-SIDE ARM]", inferred=True,
             placeholder=True), *middle,
        _seg("unpublished P7-side PCR arm", "[P7-SIDE ARM]", inferred=True,
             placeholder=True),
    ], name=f"{protocol} final library, authoritative boundary")


def _validate() -> None:
    expected = (SPACER + T7 + MARS2_TRANSCRIBED_LEADER + MARS2_RT_HANDLE
                + "N" * MARS2_WELL_NT + "N" * MARS2_UMI_NT
                + "T" * POLYT_NT + "N")
    if "".join(s.top for s in rt1_segments("MARS-seq2.0")) != expected:
        raise ValueError("MARS-seq2.0 RT1 no longer matches Supplementary Table 1")
    if POOL_READ != revcomp(POOL_ORDERED):
        raise ValueError("pool barcode read-out must be reverse-complement of the adapter")
    rt_scene("MARS-seq").rows()
    rt_scene("MARS-seq2.0").rows()


_validate()
