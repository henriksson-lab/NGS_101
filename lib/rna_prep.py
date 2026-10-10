"""Reusable molecular states for RNA-first library construction."""
from __future__ import annotations

from chemdraw import Scene, Segment, complement_segments


def seg(name: str, sequence: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name, sequence, tag, **kw)


def rna_template(name: str = "RNA insert", length: int = 34) -> list:
    """A typed RNA starting state; X is a length-preserving unknown sequence."""
    sc = Scene()
    sc.strand("RNA", [seg(name, "X" * length, placeholder=True)], label="RNA 5'→3'")
    sc.labels("RNA")
    return sc.rows()


def rna_dna_hybrid(name: str = "RNA/cDNA hybrid", length: int = 34) -> list:
    """An RNA/cDNA heteroduplex with pairing enforced by :class:`Scene`."""
    top = [seg("RNA template", "X" * length, placeholder=True)]
    return Scene.duplex(top, label=name).rows()


def single_strand(name: str, segments: list[Segment], label: str = "") -> list:
    sc = Scene()
    sc.strand(name, segments, label=label or name)
    sc.labels(name)
    return sc.rows()


def group_ii_starter_scene(rna_sequence: str, dna_sequence: str, *,
                           rna_3_block: str) -> Scene:
    """Build a blocked RNA/DNA starter with exactly one 3′ DNA overhang.

    Group-II-intron RT template switching depends on all three constraints.  Keeping
    them in the constructor prevents pages from drawing a fully paired duplex, an
    unblocked donor RNA, or a longer unspecified overhang.
    """
    if not rna_3_block.strip():
        raise ValueError("group-II starter RNA needs a named 3′ blocking group")
    if len(dna_sequence) != len(rna_sequence) + 1:
        raise ValueError("group-II starter needs exactly one overhanging DNA nucleotide")
    rna = [seg("starter RNA", rna_sequence, "r2")]
    dna = [seg("starter DNA", dna_sequence, "r2")]
    sc = Scene()
    sc.strand("starter RNA", rna, label="blocked RNA donor", mod3=rna_3_block)
    sc.anneal("starter DNA", dna, to="starter RNA",
              pair=("starter DNA", "starter RNA"), shift=-1, label="DNA primer")
    sc.labels("starter RNA")
    sc.labels("starter DNA")
    return sc
