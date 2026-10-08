"""Construct model for the individual MARS-seq protocol pages."""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, revcomp
import illumina as il
import seqprimers as sp

T7 = "TAATACGACTCACTATAGGG"
SPACER = "CGATTGAGGCCGG"
MARS2_TRANSCRIBED_LEADER = "GC"
MARS2_RT_HANDLE = "GACGTGTGCTCTTCCGATCT"
MARS2_WELL_NT = 7
MARS2_UMI_NT = 8
MARS1_WELL_NT = 6
MARS1_UMI_NT = 4
POLYT_NT = 20

POOL_ORDERED = "GACT"                 # representative ix1, as ordered
POOL_READ = revcomp(POOL_ORDERED)      # AGTC, Supplementary Table 2 "Seq barcode"
POOL_RANDOM_NT = 5
MARS1_POOL_RANDOM_NT = 3
LIG_CONSTANT = "AGATCGGAAGAGCGTCGTGTAG"

# Jaitin et al. 2014, Supplementary Tables S7--S8.  The author supplement is
# mirrored as jaitin-sm.pdf by scg_lib_structs; it is the paper's own PDF rather
# than sequence reconstructed by the upstream site.
MARS1_RT1 = (SPACER + T7 + MARS2_TRANSCRIBED_LEADER + MARS2_RT_HANDLE
             + "N" * MARS1_WELL_NT + "N" * MARS1_UMI_NT + "T" * POLYT_NT + "N")
MARS1_RT2_PRINTED = "TCTAGCCTTCTCGCAGCACATC"
MARS1_P5_PCR = il.TRUSEQ_P5_FULL
MARS1_P7_PCR = il.P7 + il.TRUSEQ_READ2

SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["R2"])


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
        return [
            _seg("5' spacer", SPACER), _seg("T7 promoter", T7, "t7"),
            _seg("transcribed leader", MARS2_TRANSCRIBED_LEADER),
            _seg("RT/PCR handle", MARS2_RT_HANDLE, "r2"),
            _seg("well barcode", "N" * MARS1_WELL_NT, "cbc", placeholder=True),
            _seg("UMI", "N" * MARS1_UMI_NT, "umi", placeholder=True),
            _seg("oligo-dT", "T" * POLYT_NT),
            _seg("N anchor", "N", placeholder=True),
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
    sc.mark("RT1", "N anchor", "poly(A) junction")
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
        adaptor = [_seg("ordered pool barcode", POOL_ORDERED, "cbc"),
                   _seg("random diversity", "N" * MARS1_POOL_RANDOM_NT,
                        placeholder=True),
                   _seg("ligation/RT2 handle", LIG_CONSTANT, "r1")]
    else:
        raise ValueError(f"unknown MARS-seq protocol: {protocol!r}")
    return Construct([
        _seg("barcoded aRNA 5' end", "[WELL BC / UMI / poly(U)]", placeholder=True),
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
        middle = [
            _seg("random diversity", "N" * MARS1_POOL_RANDOM_NT, placeholder=True),
            _seg("pool barcode", POOL_READ, "cbc"),
            _seg("sense cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
            _seg("poly(A)", "AAA", placeholder=True),
            _seg("UMI'", "N" * MARS1_UMI_NT, "umi", placeholder=True),
            _seg("well barcode'", "N" * MARS1_WELL_NT, "cbc", placeholder=True),
        ]
    else:
        raise ValueError(f"unknown MARS-seq protocol: {protocol!r}")
    if protocol == "MARS-seq":
        con = Construct([
            _seg("P5", il.P5, "p5"),
            _seg("TruSeq Read 1, non-overlap", il.TRUSEQ_READ1[4:], "r1"),
            *middle,
            _seg("TruSeq Read 2 reverse complement", revcomp(il.TRUSEQ_READ2), "r2"),
            _seg("P7 reverse complement", il.P7_RC, "p7"),
        ], name="MARS-seq final library")
        problems = sp.verify(con, SEQ_PRIMERS,
                             required_roles=("Read 1", "Read 2"))
        if problems:
            raise ValueError("invalid MARS-seq final library: " + "; ".join(problems))
        return con
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
    if "".join(s.top for s in rt1_segments("MARS-seq")) != MARS1_RT1:
        raise ValueError("MARS-seq RT1 no longer matches Supplementary Table S7")
    rt_scene("MARS-seq").rows()
    rt_scene("MARS-seq2.0").rows()


_validate()
