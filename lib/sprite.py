"""Composite identifiers and family-level split-pool tag assembly."""
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


def tag(name: str, length: int, part: int, *, group: str, role: str,
        whitelist: str) -> Segment:
    return seg(name, chr(64 + part) * length, "cbc", placeholder=True,
               feature=feature(name.lower().replace(" ", "_"), role, "combinatorial",
                               group=group, part=f"round {part}",
                               whitelist=whitelist))


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
