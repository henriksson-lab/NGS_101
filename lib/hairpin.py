"""Single-oligo hairpins and site-specific opening of their loops."""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import revcomp


@dataclass(frozen=True)
class OpenedHairpin:
    left_arm: str
    right_arm: str
    left_new_end: str
    right_new_end: str
    removed_base: str


@dataclass(frozen=True)
class Hairpin:
    name: str
    sequence: str
    stem_len: int
    overhang_3_len: int = 0

    def __post_init__(self) -> None:
        if self.stem_len < 1:
            raise ValueError("hairpin stem must contain at least one base pair")
        if self.overhang_3_len < 0:
            raise ValueError("3' overhang length cannot be negative")
        if 2 * self.stem_len + self.overhang_3_len >= len(self.sequence):
            raise ValueError("hairpin needs a non-empty loop between its stems")
        if self.left_stem != revcomp(self.right_stem):
            raise ValueError(f"{self.name}: the declared stem does not self-pair")

    @property
    def left_stem(self) -> str:
        return self.sequence[:self.stem_len]

    @property
    def right_stem(self) -> str:
        end = len(self.sequence) - self.overhang_3_len
        return self.sequence[end - self.stem_len:end]

    @property
    def loop(self) -> str:
        return self.sequence[self.stem_len:len(self.sequence)
                             - self.stem_len - self.overhang_3_len]

    @property
    def overhang_3(self) -> str:
        return self.sequence[len(self.sequence) - self.overhang_3_len:] \
            if self.overhang_3_len else ""

    def open_at(self, index: int, expected: str | None = None) -> tuple[str, str]:
        """Remove one loop nucleotide and return the two resulting linear arms."""
        if not self.stem_len <= index < len(self.sequence) - self.stem_len \
                - self.overhang_3_len:
            raise ValueError("a hairpin-opening site must lie inside the loop")
        if expected is not None and self.sequence[index] != expected:
            raise ValueError(
                f"{self.name}: expected {expected!r} at cut, found {self.sequence[index]!r}"
            )
        return self.sequence[:index], self.sequence[index + 1:]

    def user_open_at(self, index: int) -> OpenedHairpin:
        """Model USER (UDG + Endo VIII) removal of one dU from the loop."""
        left, right = self.open_at(index, expected="U")
        return OpenedHairpin(left, right, "3'-phosphate", "5'-phosphate", "dU")
