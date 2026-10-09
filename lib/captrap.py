"""Small state model for cap-trapping workflows."""
from __future__ import annotations
from dataclasses import dataclass


@dataclass(frozen=True)
class CapTrappedHybrid:
    """An RNA/cDNA hybrid eligible for streptavidin cap capture.

    Construction rejects the two common false states: a cDNA that never reached the
    cap and a cap that was not derivatized with biotin.
    """
    cdna_reaches_cap: bool
    cap_biotinylated: bool

    def __post_init__(self) -> None:
        if not self.cdna_reaches_cap:
            raise ValueError("cap trapping requires cDNA protected through the RNA cap")
        if not self.cap_biotinylated:
            raise ValueError("cap trapping requires an oxidized, biotinylated cap")

    @property
    def selected(self) -> bool:
        return True
