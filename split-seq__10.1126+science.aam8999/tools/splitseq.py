"""Construct model shared by the individual SPLiT-seq and microSPLiT pages."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from chemdraw import Construct, Scene, Segment, complement_segments

R1_HANDLE = "ACTGTGG"
R1_BARCODE = "ACTCGTAA"
R2_HEAD = "CATCGGCGTACGACT"
R2_BARCODE = "AACGTGAT"
R2_TAIL = "ATCCACGTGCTTGAG"
R3_HEAD = "CAGACGTGTGCTCTTCCGATCT"
R3_UMI_NT = 10
R3_BARCODE = "AACGTGAT"
R3_TAIL = "GTGGCCGATGTTTCG"
LINK2 = "CCACAGTCTCAAGCACG"
LINK3 = "TACGCCGATGCGAAACATCG"
TSO_BODY = "AAGCAGTGGTATCAACGCAGAGTGAAT"
TSO = TSO_BODY + "GGG"  # chemical 3' tail is rGrG+G

def _s(name, top, tag=None, **kw): return Segment(name, top, tag, **kw)

def round1(protocol: str):
    if protocol == "microSPLiT":
        return [_s("round-1 handle", R1_HANDLE), _s("barcode 1", R1_BARCODE, "cbc"),
                _s("oligo-dT", "T" * 15), _s("VN", "VN", placeholder=True)]
    if protocol == "SPLiT-seq":
        return [_s("unavailable round-1 handle", "XXXXXXXX", inferred=True, placeholder=True),
                _s("barcode 1", "XXXXXXXX", "cbc", inferred=True, placeholder=True),
                _s("oligo-dT", "T" * 15, inferred=True),
                _s("VN", "VN", inferred=True, placeholder=True)]
    raise ValueError(protocol)

def rt_scene(protocol: str):
    mrna = [_s("body", "XXXXXXXXXXXX", placeholder=True), _s("junction", "BB", placeholder=True),
            _s("poly(A)", "A" * 15)]
    sc = Scene(); sc.strand("mRNA", mrna, label="mRNA")
    sc.anneal("RT primer", round1(protocol), to="mRNA", pair=("oligo-dT", "poly(A)"),
              label="RT primer")
    sc.arrow("RT primer", "reverse transcriptase")
    return sc

def round2():
    return [_s("round-2 head", R2_HEAD), _s("barcode 2", R2_BARCODE, "cbc"),
            _s("round-2 tail", R2_TAIL)]

def round3():
    return [_s("round-3 handle", R3_HEAD, "r2"),
            _s("UMI", "N" * R3_UMI_NT, "umi", placeholder=True),
            _s("barcode 3", R3_BARCODE, "cbc"), _s("round-3 tail", R3_TAIL)]

def barcoded_first_strand(protocol: str, through: int) -> Construct:
    r1 = round1(protocol)
    segs = [*r1, _s("cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True)]
    if through >= 2:
        segs = [*round2(), *segs]
    if through >= 3:
        segs = [*round3(), *segs]
    return Construct(segs, name=f"{protocol} after barcode round {through}")

def splint_scene(round_no: int) -> Scene:
    if round_no == 2:
        left, right, source = round2(), [_s("round-1 handle", R1_HANDLE)], LINK2
        arms = [_s("round-2 tail arm", R2_TAIL[-10:]), right[0]]
    elif round_no == 3:
        left, right, source = round3(), [_s("round-2 head", R2_HEAD)], LINK3
        arms = [_s("round-3 tail arm", R3_TAIL[-10:]),
                _s("round-2 head arm", R2_HEAD[:10])]
    else: raise ValueError(round_no)
    junction = [*left, *right]
    # The source splint must equal the derived complement across the nick.
    derived = complement_segments(arms)
    if "".join(s.top for s in derived) != source:
        raise ValueError(f"round-{round_no} linker no longer bridges its two arms")
    sc = Scene(); sc.strand("joined strand", junction, label="barcoded cDNA")
    sc.anneal("linker", derived, to="joined strand",
              pair=(derived[0].name, right[0].name), label="splint")
    sc.mark("joined strand", left[-1].name, "ligated nick", through=right[0].name)
    return sc

def final_boundary(protocol: str) -> Construct:
    inner = [_s("cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
             _s("barcode chain", "[BC1 / BC2 / BC3 / UMI]", "cbc", placeholder=True,
                inferred=protocol == "SPLiT-seq")]
    return Construct([_s("unpublished left library arm", "[LEFT ARM]", inferred=True,
                         placeholder=True), *inner,
                      _s("unpublished indexed right arm", "[RIGHT ARM + i7]", inferred=True,
                         placeholder=True)], name=f"{protocol} final library boundary")

def _validate():
    rt_scene("SPLiT-seq").rows(); rt_scene("microSPLiT").rows()
    splint_scene(2).rows(); splint_scene(3).rows()
_validate()
