"""Selection states for RNA-end and RNA-modification sequencing."""
from __future__ import annotations

from batch_ngs import seg
from chemdraw import Row, Scene


def three_prime_adapter_scene(*, biotin: bool = False) -> Scene:
    parts = [seg("RNA body", "X" * 30, placeholder=True), seg("native 3-prime end", "A" * 12)]
    adapter = seg("3-prime adapter", "D" * 18, "r2", placeholder=True)
    sc = Scene(); sc.strand("RNA", [*parts, adapter], label="3-prime-adapter-ligated RNA")
    sc.junction("RNA", "native 3-prime end", "3-prime adapter", "RNA ligation")
    if biotin:
        sc.mark("RNA", "3-prime adapter", "biotin capture handle")
    sc.labels("RNA")
    return sc


def terminal_fragment_selection_rows() -> list[Row]:
    return [
        Row(chunks=[("adapter-ligated RNA — partial RNase T1 digest", None, False)]),
        Row(chunks=[("streptavidin retains only fragments carrying the native 3-prime end", "w1", False)]),
    ]


def immunoprecipitation_rows(mark: str = "m6A") -> list[Row]:
    return [
        Row(chunks=[(f"fragmented poly(A) RNA ---- ● {mark} ---- RNA", None, False)]),
        Row(chunks=[(f"                         ^ anti-{mark} antibody", "cbc", False)]),
        Row(chunks=[("                         └─ enriched RNA fragment", "w1", False)]),
    ]
