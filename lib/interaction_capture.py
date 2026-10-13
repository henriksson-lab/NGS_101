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


@dataclass(frozen=True)
class FixedRNPContact:
    """Two RNA fragments retained in one cross-linked RNP neighbourhood.

    Unlike :class:`ReversibleDuplexCrosslink`, this model deliberately does not anneal
    the two RNA strands: RIC-seq proximity can be protein-mediated and therefore does
    not establish Watson-Crick pairing.  Matching contact labels show that both strands
    belong to one fixed neighbourhood without inventing a base-resolved fold.
    """

    crosslinker: str

    def __post_init__(self) -> None:
        if not self.crosslinker.strip():
            raise ValueError("an RNP contact needs a named crosslinking mechanism")

    def complex(self) -> Scene:
        sc = Scene()
        sc.strand("RNA A", [_part("RNA A fragment", "X"),
                            _part("RNP contact A", "P", 8)],
                  label="RNA fragment A")
        sc.strand("RNA B", [_part("RNP contact B", "Q", 8),
                            _part("RNA B fragment", "Y")],
                  label="RNA fragment B")
        label = f"{self.crosslinker}: shared fixed RNP neighbourhood"
        sc.mark("RNA A", "RNP contact A", label)
        sc.mark("RNA B", "RNP contact B", label)
        sc.note("RNA B", "proximity is asserted; RNA–RNA base pairing is not")
        return sc

    def chimera(self) -> Scene:
        parts = [_part("RNA A fragment", "X"), _part("RNA B fragment", "Y")]
        sc = Scene(); sc.strand("chimera", parts, label="RNA proximity-ligation chimera")
        sc.junction("chimera", parts[0].name, parts[1].name, "RNA proximity ligation")
        sc.labels("chimera")
        return sc


@dataclass(frozen=True)
class BivalentRNAChromatinBridge:
    """Role-level two-ended linker connecting RNA to nearby genomic DNA.

    The helper separates the first RNA-end ligation from the second DNA-end ligation,
    preventing a page from depicting the bridge as an unexplained linear word between
    the two analytes.  Role tokens are used because linker bases are protocol-specific.
    """

    name: str

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise ValueError("a bivalent RNA–DNA bridge needs a name")

    def unjoined(self) -> Scene:
        sc = Scene()
        sc.strand("RNA", [_part("chromatin-associated RNA", "R")], label="nearby RNA")
        sc.strand("bridge", [_part("RNA-reactive linker end", "L", 10, "me"),
                             _part("annealed linker body", "M", 14, "me"),
                             _part("DNA-reactive linker end", "N", 10, "me")],
                  label=self.name)
        sc.strand("DNA", [_part("nearby genomic DNA", "D")], label="nearby chromatin")
        sc.mark("bridge", "RNA-reactive linker end", "first ligation: RNA end")
        sc.mark("bridge", "DNA-reactive linker end", "second ligation: DNA end")
        sc.note("DNA", "crosslinked nuclear proximity; no RNA–DNA base pairing asserted")
        return sc

    def rna_joined(self) -> Scene:
        joined = [_part("chromatin-associated RNA", "R"),
                  _part("RNA-reactive linker end", "L", 10, "me"),
                  _part("annealed linker body", "M", 14, "me"),
                  _part("free DNA-reactive linker end", "N", 10, "me")]
        sc = Scene(); sc.strand("RNA–linker", joined, label="RNA-end-ligated bridge")
        sc.junction("RNA–linker", joined[0].name, joined[1].name, "RNA-end ligation")
        sc.strand("DNA", [_part("nearby genomic DNA", "D")], label="nearby chromatin")
        sc.mark("RNA–linker", "free DNA-reactive linker end", "available for second ligation")
        return sc

    def chimera(self) -> Scene:
        parts = [_part("chromatin-associated RNA", "R"),
                 _part(self.name, "L", 24, "me"),
                 _part("nearby genomic DNA", "D")]
        sc = Scene(); sc.strand("contact", parts, label="RNA–linker–DNA contact chimera")
        sc.junction("contact", parts[0].name, parts[1].name, "RNA-end ligation")
        sc.junction("contact", parts[1].name, parts[2].name, "DNA-end ligation")
        sc.labels("contact")
        return sc
