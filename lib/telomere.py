"""Shared structures for chromosome-end tagging and telomere capture workflows."""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import Scene, Segment, complement_segments, revcomp

HUMAN_G_REPEAT = "TTAGGG"
HUMAN_C_REPEAT = revcomp(HUMAN_G_REPEAT)


def phases(sequence: str) -> tuple[str, ...]:
    """Every cyclic phase of a repeat, in deterministic order."""
    sequence = sequence.upper()
    if not sequence or set(sequence) - set("ACGT"):
        raise ValueError("repeat must be non-empty DNA")
    return tuple(sequence[i:] + sequence[:i] for i in range(len(sequence)))


HUMAN_G_PHASES = phases(HUMAN_G_REPEAT)
HUMAN_C_PHASES = phases(HUMAN_C_REPEAT)


@dataclass(frozen=True)
class TelomereEnd:
    """A duplex telomere ending in a native G-rich 3-prime overhang."""

    repeat: str = HUMAN_G_REPEAT
    duplex_copies: int = 3
    overhang_copies: int = 2
    subtelomere_bases: int = 16

    def __post_init__(self) -> None:
        if self.duplex_copies < 1 or self.overhang_copies < 1:
            raise ValueError("telomere drawing needs duplex repeat and overhang")
        if self.subtelomere_bases < 1:
            raise ValueError("telomere drawing needs a subtelomeric anchor")

    @property
    def duplex(self) -> str:
        return self.repeat * self.duplex_copies

    @property
    def overhang(self) -> str:
        return self.repeat * self.overhang_copies

    def scene(self) -> Scene:
        top = [Segment("subtelomere", "X" * self.subtelomere_bases, placeholder=True),
               Segment("duplex telomere", self.duplex, "r1"),
               Segment("G-rich 3-prime overhang", self.overhang, "r2")]
        bottom = complement_segments(top[:2])
        sc = Scene(); sc.strand("G-rich strand", top, label="chromosome")
        sc.anneal("C-rich strand", bottom, to="G-rich strand",
                  pair=("duplex telomere'", "duplex telomere"), label="chromosome")
        sc.mark("G-rich strand", "G-rich 3-prime overhang", "native terminal overhang")
        return sc


@dataclass(frozen=True)
class PhasedTelomereAdapter:
    """Adapter family whose terminal arms cover every repeat phase."""

    core: str
    arm_repeat: str = HUMAN_C_REPEAT
    arm_copies: int = 3
    phosphorylated_5: bool = True

    def __post_init__(self) -> None:
        if not self.core or set(self.core.upper()) - set("ACGTNX"):
            raise ValueError("adapter core must be DNA or explicit N/X placeholders")
        if not self.arm_repeat or set(self.arm_repeat.upper()) - set("ACGT"):
            raise ValueError("telomere arm repeat must be real DNA")
        if self.arm_copies < 1:
            raise ValueError("adapter needs at least one repeat copy")

    @property
    def arms(self) -> tuple[str, ...]:
        return tuple((p + self.arm_repeat * self.arm_copies)[:
                     len(self.arm_repeat) * self.arm_copies]
                     for p in phases(self.arm_repeat))

    def annealed_scene(self, end: TelomereEnd | None = None, phase: int = 0) -> Scene:
        end = end or TelomereEnd()
        arm = self.arms[phase]
        top = [Segment("subtelomere", "X" * end.subtelomere_bases, placeholder=True),
               Segment("duplex telomere", end.duplex, "r1"),
               Segment("G-rich overhang", end.overhang, "r2")]
        bottom = complement_segments(top[:2])
        sc = Scene(); sc.strand("G-rich strand", top, label="chromosome")
        sc.anneal("C-rich strand", bottom, to="G-rich strand",
                  pair=("duplex telomere'", "duplex telomere"), label="chromosome")
        adapter = [Segment("adapter core", self.core, "r3", placeholder="N" in self.core or "X" in self.core),
                   Segment("telomere-complement arm", arm, "r2")]
        sc.anneal("terminal adapter", adapter, to="G-rich strand",
                  pair=("telomere-complement arm", "G-rich overhang"),
                  shift=len(end.overhang) - len(arm), label="adapter",
                  unpaired=("adapter core",), mod5="p" if self.phosphorylated_5 else "")
        sc.footer("** ligase seals the adapter 5′ end to the native C-rich 3′ end",
                  "terminal adapter", "adapter core")
        return sc


@dataclass(frozen=True)
class AffinityCapture:
    """A tagged molecule selected on a solid phase and released by a named cut."""

    affinity: str
    release_enzyme: str | None = None

    def steps(self) -> tuple[str, ...]:
        out = (f"bind {self.affinity}-tagged molecules", "wash away untagged DNA")
        return out + ((f"release with {self.release_enzyme}",) if self.release_enzyme else ())
