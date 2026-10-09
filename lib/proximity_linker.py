"""Double-ended bridge linkers for proximity ligation."""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import Scene, Segment, revcomp


@dataclass(frozen=True)
class ThreePrimeBridgeLinker:
    """Duplex bridge with one 3-prime overhang at each end.

    Both oligos are supplied 5' to 3'.  Their terminal ``overhang`` bases are removed
    before duplex compatibility is checked, making the advertised two-ended topology a
    construction invariant rather than a caption.
    """

    name: str
    top_5to3: str
    bottom_5to3: str
    overhang: str = "T"
    modified_top_position: int | None = None

    def __post_init__(self) -> None:
        if not self.overhang or not self.top_5to3.endswith(self.overhang):
            raise ValueError("top oligo does not end in the declared 3-prime overhang")
        if not self.bottom_5to3.endswith(self.overhang):
            raise ValueError("bottom oligo does not end in the declared 3-prime overhang")
        n = len(self.overhang)
        if self.top_5to3[:-n] != revcomp(self.bottom_5to3[:-n]):
            raise ValueError("bridge-linker duplex cores are not reverse complements")
        if self.modified_top_position is not None:
            if not 0 <= self.modified_top_position < len(self.top_5to3) - n:
                raise ValueError("the internal modification must lie in the duplex core")

    @property
    def top_core(self) -> str:
        return self.top_5to3[:-len(self.overhang)]

    @property
    def bottom_core(self) -> str:
        return self.bottom_5to3[:-len(self.overhang)]

    def scene(self) -> Scene:
        split = self.modified_top_position
        if split is None:
            top_core = [Segment("bridge core", self.top_core, "cbc")]
        else:
            top_core = [Segment("bridge core left", self.top_core[:split], "cbc"),
                        Segment("modified bridge base", self.top_core[split], "umi"),
                        Segment("bridge core right", self.top_core[split + 1:], "cbc")]
        top = [*top_core, Segment("top 3-prime overhang", self.overhang, "r2")]
        bottom = [Segment("bottom bridge core", self.bottom_core, "cbc"),
                  Segment("bottom 3-prime overhang", self.overhang, "r2")]
        sc = Scene()
        sc.strand("bridge top", top, label=self.name, mod5="p")
        sc.anneal("bridge bottom", bottom, to="bridge top",
                  pair=("bottom bridge core", top_core[0].name),
                  unpaired=("top 3-prime overhang", "bottom 3-prime overhang"),
                  label=self.name, mod5="p")
        if split is not None:
            sc.mark("bridge top", "modified bridge base", "internal biotin")
        return sc
