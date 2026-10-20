"""Reusable molecular states for chromatin and epigenetic library protocols."""
from __future__ import annotations

from batch_ngs import seg
from chemdraw import Row, Scene


def pbat_priming_scene(round_name: str, tail_name: str, tail_seq: str,
                       *, mod5: str = "") -> Scene:
    template = [seg("bisulfite-converted template", "X" * 26, placeholder=True)]
    primer = [seg(tail_name, tail_seq, "r1"),
              seg("random N4", "N" * 4, placeholder=True)]
    sc = Scene(); sc.strand("template", template, label=round_name + " template")
    sc.anneal("primer", primer, to="template",
              pair=("random N4", "bisulfite-converted template"),
              shift=len(template[0]) - len(primer[-1]),
              unpaired=("random N4",), label=round_name + " random primer", mod5=mod5)
    sc.mark("primer", tail_name, "5′ sequencing-adapter tail")
    sc.arrow("primer", "polymerase extension")
    return sc


def antibody_enzyme_rows(enzyme: str, action: str) -> list[Row]:
    return [
        Row(chunks=[("chromatin ---- [target protein or histone mark] ---- chromatin", None, False)]),
        Row(chunks=[("                         ^ primary antibody", "cbc", False)]),
        Row(chunks=[(f"                         | Protein A–{enzyme}", "me", False)]),
        Row(chunks=[(f"                         └─ {action}", None, False)]),
    ]


def dimelo_mark_rows() -> list[Row]:
    return antibody_enzyme_rows("Hia5", "SAM-dependent m6A marks deposited on nearby DNA") + [
        Row(chunks=[("single long molecule: ---A--m6A---A----m6A------A---", "w1", False)]),
    ]


def strand_selection_rows() -> list[Row]:
    return [
        Row(chunks=[("parental template       5′ -------- parental -------- 3′", None, False)]),
        Row(chunks=[("new BrdU strand         3′ --BrdU--BrdU--BrdU------ 5′", "w1", False)]),
        Row(chunks=[("                                      UV + Hoechst ↓", None, False)]),
        Row(chunks=[("retained library strand 5′ -------- parental -------- 3′", None, False)]),
    ]
