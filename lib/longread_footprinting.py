"""Shared molecular states for adenine-methyltransferase fiber assays.

The accessible DNA is marked before extraction.  Library construction must remain
PCR-free because copying would replace the modified template with unmodified DNA.
"""
from __future__ import annotations

from chemdraw import Construct, Row, Scene, Segment
from dumbbell import Dumbbell, dumbbell_rows


def marked_fiber(enzyme: str, *, target: str) -> list[Row]:
    """A duplex fiber with methyltransferase access represented as anchored marks."""
    con = Construct([
        Segment("protein-protected DNA", "X" * 14, placeholder=True),
        Segment(target, "A" * 12, "w1"),
        Segment("protein-protected DNA 2", "X" * 14, placeholder=True),
        Segment(target + " 2", "A" * 12, "w1"),
        Segment("protein-protected DNA 3", "X" * 14, placeholder=True),
    ], name="chromatin fiber")
    sc = Scene.duplex(list(con), label="native chromatin DNA")
    sc.mark("top", target, f"{enzyme}: accessible A -> m6dA")
    sc.mark("top", target + " 2", f"{enzyme}: accessible A -> m6dA")
    sc.note("top", "protein-bound intervals remain comparatively unmarked")
    return sc.rows()


def conventional_smrtbell(insert_name: str | Construct) -> Dumbbell:
    """A repaired/dA-tailed insert closed by two proprietary PacBio hairpins."""
    insert = (insert_name if isinstance(insert_name, Construct) else
              Construct([Segment(insert_name, "X" * 54, placeholder=True)], name=insert_name))
    return Dumbbell(insert, "N" * 18, "N" * 18, name=insert.name + " SMRTbell",
                    insert_overhang="A", adapter_overhang="T")


def hairpin_tagmented_smrtbell(
        insert_name: str, *,
        adapter_sequence: str | None = None, stem_length: int = 0) -> Dumbbell:
    """The final closed product of hairpin-Tn5 transfer followed by gap repair.

    Both hairpins are mandatory in :class:`Dumbbell`; this prevents a caller from
    presenting a one-ended tagmentation product as a sequenceable PacBio molecule.
    """
    insert = Construct([Segment(insert_name, "X" * 54, placeholder=True)], name=insert_name)
    if adapter_sequence is None:
        return Dumbbell(insert, "N" * 18, "N" * 18,
                        name="hairpin-tagmented SMRTbell",
                        insert_overhang="A", adapter_overhang="T")
    return Dumbbell(insert, adapter_sequence, adapter_sequence,
                    name="hairpin-tagmented SMRTbell", insert_overhang="",
                    adapter_overhang="", left_stem_length=stem_length,
                    right_stem_length=stem_length)


def smrtbell_rows(molecule: Dumbbell) -> list[Row]:
    return dumbbell_rows(molecule)


def pacbio_entry(molecule: Dumbbell, *, adapter_sequence_known: bool = False) -> str:
    """Source-safe description of primer/polymerase geometry for a closed template."""
    adapter = ("either published hairpin adapter" if adapter_sequence_known
               else "either proprietary hairpin adapter")
    return (f"A proprietary PacBio sequencing primer anneals in {adapter}. "
            "The bound polymerase traverses the marked insert, the opposite hairpin, "
            "the reverse-complement insert and the starting hairpin repeatedly to form "
            "a circular-consensus read.")
