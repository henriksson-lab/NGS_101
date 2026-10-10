"""Checked molecular states for transcriptome-wide RNA interaction assays.

The important distinction is made in the constructor: psoralen assays retain a
base-paired RNA duplex through enrichment, whereas protein-mediated assays retain two
RNAs on a cross-linked RNP.  Both become one covalent RNA chimera only at the explicit
proximity-ligation transition.
"""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import Scene, Segment


def _part(name: str, base: str = "X", n: int = 24, tag: str | None = None) -> Segment:
    return Segment(name, base * n, tag, placeholder=True)


@dataclass(frozen=True)
class ReversibleDuplexCrosslink:
    reagent: str
    selection: str
    reversible: bool = True

    def __post_init__(self) -> None:
        if not self.reagent.strip() or not self.selection.strip():
            raise ValueError("RNA-duplex crosslinking needs a reagent and selection route")
        if not self.reversible:
            raise ValueError("proximity ligation requires a reversible duplex crosslink")

    def duplex(self) -> Scene:
        top = [_part("RNA partner A", "X")]
        bottom = [_part("RNA partner B", "x")]
        sc = Scene(); sc.strand("RNA A", top, label="base-paired RNA")
        sc.anneal("RNA B", bottom, to="RNA A", pair=("RNA partner B", "RNA partner A"),
                  label="base-paired RNA")
        sc.mark("RNA A", "RNA partner A", f"{self.reagent} interstrand crosslink")
        sc.note("RNA B", self.selection)
        return sc

    def chimera(self) -> Scene:
        parts = [_part("RNA partner A fragment", "X"),
                 _part("RNA partner B fragment", "Y")]
        sc = Scene(); sc.strand("chimera", parts, label="RNA proximity-ligation product")
        sc.junction("chimera", parts[0].name, parts[1].name, "RNA proximity ligation")
        sc.note("chimera", f"reverse {self.reagent} crosslink after ligation")
        sc.labels("chimera")
        return sc


@dataclass(frozen=True)
class ProteinMediatedRNAContact:
    linker_name: str
    selectable: bool

    def __post_init__(self) -> None:
        if not self.linker_name.strip() or not self.selectable:
            raise ValueError("protein-mediated RNA capture needs a selectable linker")

    def complex(self) -> Scene:
        sc = Scene(); sc.strand("RNA A", [_part("protein-bound RNA A")], label="cross-linked RNP")
        sc.strand("RNA B", [_part("protein-bound RNA B", "Y")], label="cross-linked RNP")
        sc.note("RNA A", "both RNAs are retained on one UV-cross-linked protein complex")
        return sc

    def chimera(self) -> Scene:
        parts = [_part("RNA A fragment"), _part(self.linker_name, "L", 16, "cbc"),
                 _part("RNA B fragment", "Y")]
        sc = Scene(); sc.strand("chimera", parts, label="selectable RNA–linker–RNA chimera")
        sc.junction("chimera", parts[0].name, parts[1].name, "RNA ligation")
        sc.junction("chimera", parts[1].name, parts[2].name, "RNA ligation")
        sc.mark("chimera", self.linker_name, "affinity selection")
        sc.labels("chimera")
        return sc
