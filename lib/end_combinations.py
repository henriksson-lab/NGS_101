"""Two-ended product classes derived from a pool of loaded adapter ends.

This is an end-identity model. It does not model transposition efficiency, fragment
frequency, or whether a PCR with a single primer species can make any product.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations_with_replacement

from pcr import PCRPrimer


@dataclass(frozen=True)
class LoadedEnd:
    name: str
    sequence: str  # 5' entry site through the shared transposase-binding end

    def __post_init__(self) -> None:
        if not self.name.strip() or not self.sequence.strip():
            raise ValueError("loaded end needs a name and sequence")


@dataclass(frozen=True)
class LoadedEndPool:
    ends: tuple[LoadedEnd, ...]
    shared_anchor: str

    def __post_init__(self) -> None:
        if len(self.ends) < 2 or len({end.name for end in self.ends}) != len(self.ends):
            raise ValueError("end pool needs at least two differently named ends")
        if not self.shared_anchor or any(
                not end.sequence.endswith(self.shared_anchor) or
                len(end.sequence) <= len(self.shared_anchor) for end in self.ends):
            raise ValueError("each loaded end needs an entry site before the shared anchor")

    def pairs(self) -> tuple[tuple[str, str], ...]:
        return tuple(combinations_with_replacement(
            sorted(end.name for end in self.ends), 2))


@dataclass(frozen=True)
class EndProductGroup:
    output: str
    pairs: tuple[tuple[str, str], ...]

    def __post_init__(self) -> None:
        if not self.output or not self.pairs:
            raise ValueError("end-product group needs an output and end pairs")


@dataclass(frozen=True)
class EndCombinationGeometry:
    primers: tuple[PCRPrimer, PCRPrimer]
    groups: tuple[EndProductGroup, ...]
    pcr_step: str

    def __post_init__(self) -> None:
        if len({primer.sequence for primer in self.primers}) != 2:
            raise ValueError("two distinct primer species are needed")
        if not self.groups:
            raise ValueError("end combinations need declared product groups")
        if not self.pcr_step:
            raise ValueError("end combinations need the downstream PCR step")


def validate_groups(pool: LoadedEndPool, geometry: EndCombinationGeometry,
                    output_is_waste: dict[str, bool]) -> str | None:
    """Return a mismatch, or None when all possible end classes are accounted for.

    A retained class needs one distinct primer's exact 3' footprint at each end.
    Rejected classes are those missing that pair of distinct entry sites. The
    classification concerns eligibility for this two-primer library, not absolute
    amplification or product abundance.
    """
    names = {end.name for end in pool.ends}
    sites = {}
    for end in pool.ends:
        entry = end.sequence[:-len(pool.shared_anchor)]
        hits = tuple(i for i, primer in enumerate(geometry.primers)
                     if entry.startswith(primer.sequence[-primer.anneal_nt:]))
        if len(hits) > 1:
            return f"loaded end {end.name!r} matches more than one primer footprint"
        sites[end.name] = hits[0] if hits else None

    expected = set(pool.pairs())
    observed = set()
    for group in geometry.groups:
        if group.output not in output_is_waste:
            return f"end-product output {group.output!r} is absent from the step"
        for pair in group.pairs:
            if len(pair) != 2 or any(name not in names for name in pair):
                return f"end pair {pair!r} names an unknown loaded end"
            canonical = tuple(sorted(pair))
            if canonical in observed:
                return f"end pair {canonical!r} is assigned more than once"
            observed.add(canonical)
            hits = (sites[canonical[0]], sites[canonical[1]])
            eligible = set(hits) == {0, 1}
            is_waste = output_is_waste[group.output]
            if eligible == is_waste:
                return (f"end pair {canonical!r} is assigned to {group.output!r}, "
                        f"but {'has' if eligible else 'lacks'} both distinct primer sites")
    if observed != expected:
        return f"end classes are missing: {sorted(expected - observed)!r}"
    return None
