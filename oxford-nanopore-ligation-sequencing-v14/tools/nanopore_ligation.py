"""Molecular model for Oxford Nanopore SQK-LSK114 ligation sequencing."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment
from endprep import (ThreePrimeOverhangAdapter, dA_tailed_scene, ligate_both_ends,
                     repair_and_dA_tail)


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


LA_TOP = "CCTGTACTTCGTTCAGTTACGTATTGCT"
LA_BOTTOM = "GCAATACGTAACTGAACGAAGTACAGG"
LIGATION_ADAPTER = ThreePrimeOverhangAdapter(
    "SQK-LSK114 Ligation Adapter disclosed duplex", LA_TOP, LA_BOTTOM
)
INSERT = Construct([seg("genomic DNA", "X" * 42, placeholder=True)], name="DNA fragment")
END_PREP = repair_and_dA_tail(INSERT)
SEALED_JUNCTIONS = ligate_both_ends(END_PREP, LIGATION_ADAPTER)
END_PREP_TIMES = ((20, 5), (65, 5))


def end_prep_scene() -> Scene:
    return dA_tailed_scene(END_PREP)


def adapter_scene() -> Scene:
    top = [seg("27-bp adapter core", LIGATION_ADAPTER.core, "r1"),
           seg("3-prime dT", LIGATION_ADAPTER.overhang_3, "r2")]
    bottom = [seg("adapter bottom", LIGATION_ADAPTER.bottom, "r1")]
    sc = Scene(); sc.strand("adapter top", top, label="LA top")
    sc.anneal("adapter bottom", bottom, to="adapter top",
              pair=("adapter bottom", "27-bp adapter core"), label="LA bottom")
    sc.mark("adapter top", "3-prime dT", "pairs the DNA dA")
    return sc


def final_library() -> Construct:
    left = SEALED_JUNCTIONS.left
    right = SEALED_JUNCTIONS.right
    return Construct([
        seg("left motor-loaded leader", "[motor-loaded leader]", "r2", placeholder=True),
        seg("left adapter core", LIGATION_ADAPTER.core, "r1"),
        seg("left T:A junction", left.adapter_overhang, bottom=left.insert_overhang),
        seg("genomic DNA", "X" * 42, placeholder=True),
        seg("right A:T junction", right.insert_overhang, bottom=right.adapter_overhang),
        seg("right adapter core", LIGATION_ADAPTER.bottom, "r1"),
        seg("right motor-loaded leader", "[motor-loaded leader]", "r2", placeholder=True),
    ], name="SQK-LSK114 adapted duplex")


def final_scene() -> Scene:
    lib = final_library()
    sc = Scene.duplex(list(lib), label="adapted DNA",
                      unpaired=("left motor-loaded leader'", "right motor-loaded leader'"))
    sc.junction("top", "left T:A junction", "genomic DNA", "sealed ligation")
    sc.junction("top", "genomic DNA", "right A:T junction", "sealed ligation")
    return sc
