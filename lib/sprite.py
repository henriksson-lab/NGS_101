"""Split-pool recognition-by-tag-extension constructs.

SPRITE identifiers are composite: in scSPRITE the first three tag parts identify a
nucleus while all six identify a spatial complex.  ``CompositeIdentifier`` records that
overlap without forcing one DNA segment to pretend it has only one biological meaning.
"""
from __future__ import annotations

from dataclasses import dataclass

from batch_ngs import joined_scene, seg, truseq_library
from chemdraw import Construct, Row, Segment, feature


@dataclass(frozen=True)
class CompositeIdentifier:
    id: str
    role: str
    segments: tuple[str, ...]

    def bind(self, construct: Construct) -> "CompositeIdentifier":
        names = {s.name for s in construct}
        missing = [n for n in self.segments if n not in names]
        if missing:
            raise ValueError(f"{self.id}: missing barcode parts {missing}")
        if len(set(self.segments)) != len(self.segments):
            raise ValueError(f"{self.id}: a barcode part is repeated")
        return self


def tag(name: str, length: int, part: int, *, group: str, role: str) -> Segment:
    return seg(name, chr(64 + part) * length, "cbc", placeholder=True,
               feature=feature(name.lower().replace(" ", "_"), role, "combinatorial",
                               group=group, part=f"round {part}",
                               whitelist="published 96-well SPRITE tag set"))


def bulk_barcode() -> tuple[list[Segment], CompositeIdentifier]:
    parts = [tag("DPM tag", 9, 1, group="sprite_cluster", role="inline_barcode"),
             tag("Odd tag 1", 17, 2, group="sprite_cluster", role="inline_barcode"),
             tag("Even tag", 17, 3, group="sprite_cluster", role="inline_barcode"),
             tag("Odd tag 2", 17, 4, group="sprite_cluster", role="inline_barcode"),
             tag("Terminal tag", 9, 5, group="sprite_cluster", role="inline_barcode")]
    return parts, CompositeIdentifier("sprite_cluster", "spatial complex", tuple(x.name for x in parts))


def single_cell_barcode() -> tuple[list[Segment], tuple[CompositeIdentifier, ...]]:
    names = ("DPM tag", "Odd cell tag", "Even cell tag", "Odd spatial tag",
             "Even spatial tag", "Y-even terminal tag")
    lengths = (9, 17, 17, 17, 17, 9)
    parts=[]
    for i,(name,n) in enumerate(zip(names,lengths),1):
        role="cell_barcode" if i <= 3 else "inline_barcode"
        group="cell_id" if i <= 3 else "spatial_extension"
        parts.append(tag(name,n,i,group=group,role=role))
    ids=(CompositeIdentifier("cell_id","cell",tuple(x.name for x in parts[:3])),
         CompositeIdentifier("spatial_cluster","spatial complex",tuple(x.name for x in parts)))
    return parts,ids


def ligation_round_rows(parts: list[Segment]) -> list[Row]:
    chain = [seg("genomic fragment", "X" * 30, placeholder=True), *parts]
    boundaries = tuple((chain[i].name, chain[i + 1].name, f"ligation round {i + 1}")
                       for i in range(len(chain) - 1))
    return joined_scene(chain, boundaries, label="split-pool tag chain", duplex=False).rows()


def final_library(parts: list[Segment], name: str):
    """Place tag parts so R1 sees DPM then genome and R2 sees remaining tags in order."""
    insert=[parts[0],seg("genomic DNA", "X"*42, placeholder=True),*reversed(parts[1:])]
    lib,primers=truseq_library(insert,name)
    return lib,primers
