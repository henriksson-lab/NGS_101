"""
Padlock / molecular inversion probes: capture, gap-fill, circularisation.

A padlock probe is one linear ssDNA whose two ends are complementary to two places on a
target, a short distance apart. Both ends anneal, the gap between them is filled by a
polymerase, and a ligase joins the new 3' end to the probe's 5' phosphate. The result is a
covalently closed circle that contains both the probe's own backbone and a copy of the
captured target.

Three properties follow, and they are the whole reason the chemistry is used:

* **The circle is the selection.** Probes that failed to find a target, or failed to
  ligate, stay linear -- and a 3'->5' exonuclease destroys everything linear. No gel, no
  bead cleanup: the enzymology throws away the failures.
* **It is an inversion.** The captured target ends up flanked by probe backbone, so the
  probe's UMI and index travel with the captured molecule.
* **Counting, not amplifying.** A UMI in the backbone is attached before any amplification,
  so the readout counts original molecules instead of PCR products.

A "molecular inversion probe" (MIP) is a padlock probe used this way for quantification.
"""

from __future__ import annotations

from dataclasses import dataclass

from chemdraw import revcomp, tm
from plasmid import Plasmid, find_both


@dataclass
class Padlock:
    """A padlock/MIP probe, written 5'->3' as:  lig_arm - backbone - ext_arm.

    ext_arm   the probe's 3' end. Anneals UPSTREAM of the gap; the polymerase extends
              from it, across the gap, toward the ligation arm.
    lig_arm   the probe's 5' end, carrying the phosphate. Anneals DOWNSTREAM of the gap
              and receives the ligation that closes the circle.
    backbone  not complementary to the target: UMI, primer sites, index.

    Arm sequences are given in the same sense as the target strand they are drawn
    against, which is how probe tables are normally published.
    """
    name: str
    ext_arm: str
    lig_arm: str
    backbone: str
    phosphorylated: bool = True

    def sequence(self) -> str:
        return self.lig_arm + self.backbone + self.ext_arm

    def __len__(self) -> int:
        return len(self.sequence())

    def arm_tms(self) -> tuple[float, float]:
        return tm(self.ext_arm), tm(self.lig_arm)

    def arm_tm_gap(self) -> float:
        """Ligation arm minus extension arm.

        Published designs make the ligation arm the hotter of the two, so that it stays
        bound while the polymerase approaches it -- otherwise the polymerase displaces it
        and the circle never closes.
        """
        e, l = self.arm_tms()
        return l - e


@dataclass
class Capture:
    """One successful capture event."""
    probe: Padlock
    target: str
    ext_end: int          # 0-based, first base after the extension arm
    lig_start: int        # 0-based start of the ligation arm
    gap: int              # bases the polymerase must synthesise
    fill: str             # the captured target sequence
    circle: Plasmid       # the closed product

    @property
    def captured_span(self) -> int:
        """Target bases spanned, arms included."""
        return len(self.probe.ext_arm) + self.gap + len(self.probe.lig_arm)


def capture(target: Plasmid, probe: Padlock, max_gap: int = 2000) -> list[Capture]:
    """Find where this probe closes on this target.

    Both arms must land on the same strand, in the order ext -> lig, with a gap the
    polymerase can plausibly cross. Circular targets are handled, so a capture spanning
    the origin is found like any other.
    """
    out: list[Capture] = []
    n = len(target)
    for e_pos, e_strand in find_both(target.seq, probe.ext_arm, target.circular):
        for l_pos, l_strand in find_both(target.seq, probe.lig_arm, target.circular):
            if e_strand != l_strand:
                continue
            if e_strand == 1:
                e_end, l_start = (e_pos + len(probe.ext_arm)) % n, l_pos
            else:
                # on the minus strand the probe runs the other way along the plus-strand map
                e_end, l_start = (l_pos + len(probe.lig_arm)) % n, e_pos
            gap = (l_start - e_end) % n
            if not 0 <= gap <= max_gap:
                continue
            fill = target.sub(e_end, e_end + gap)
            if e_strand == -1:
                fill = revcomp(fill)
            # the circle: captured target, then the probe backbone, closed on itself
            seq = probe.ext_arm + fill + probe.lig_arm + probe.backbone
            out.append(Capture(probe, target.name, e_end, l_start, gap, fill,
                               Plasmid(f"{probe.name}o{target.name}", seq, True, [])))
    return sorted(out, key=lambda c: c.gap)


def survives_exonuclease(molecule: Plasmid | str) -> bool:
    """Exonuclease I/III digest everything with a free end. Only closed circles remain.

    This is the step that makes the method specific: an unreacted probe, a probe that
    annealed but never ligated, and the genomic DNA itself are all linear.
    """
    return isinstance(molecule, Plasmid) and molecule.circular
