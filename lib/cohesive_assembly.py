"""Geometry of explicitly supplied DNA fragments and cohesive-end polarities.

This module chooses no enzyme, biological part, fusion vocabulary or assembly
order. It checks one declared circular product, not possible-product frequencies,
reaction efficiency or pool coverage. Release from source plasmids is a separate
relation: a fragment declaration alone does not establish a restriction digest.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

from chemdraw import Construct, Row, Scene, Segment, complement_segments, revcomp


@dataclass(frozen=True)
class CohesiveFragment:
    """A double-stranded core between two source-declared cohesive ends.

    Both fusion strings use the same reference top-strand axis. The physical
    left 5-prime overhang is on top, while a left 3-prime overhang is on the
    bottom strand. At the right end these assignments reverse. Bottom-strand
    overhangs are written as reverse complements of the supplied reference axis.
    """

    name: str
    left: str
    core: tuple[Segment, ...]
    right: str
    left_polarity: str
    right_polarity: str

    def __post_init__(self) -> None:
        if not self.name.strip() or not self.core:
            raise ValueError("a cohesive fragment needs a name and double-stranded core")
        for end in (self.left, self.right):
            if not end or set(end) - set("ACGT"):
                raise ValueError("cohesive ends must have explicit, uppercase DNA bases")
        if any(p not in {"5-prime", "3-prime"} for p in (self.left_polarity, self.right_polarity)):
            raise ValueError("each cohesive end needs an explicit 5-prime or 3-prime polarity")
        Construct(self.core, self.name)  # reject conflicting segment/feature identities
        if not all(seg.name and seg.top for seg in self.core):
            raise ValueError("cohesive fragment core segments need names and content")
        self.scene()  # validate any caller-supplied bottom strand at construction

    def scene(self) -> Scene:
        core = tuple(replace(seg, name="core " + seg.name) for seg in self.core)
        bottom = complement_segments(core)
        sc = Scene()
        top_parts, bottom_parts = list(core), list(bottom)
        if self.left_polarity == "5-prime":
            top_parts.insert(0, Segment("left fusion", self.left))
        else:
            bottom_parts.append(Segment("left fusion bottom", revcomp(self.left)))
        if self.right_polarity == "5-prime":
            bottom_parts.insert(0, Segment("right fusion bottom", revcomp(self.right)))
        else:
            top_parts.append(Segment("right fusion", self.right))
        sc.strand("top", top_parts, label=self.name)
        sc.anneal("bottom", bottom_parts,
                  to="top", pair=(bottom[-1].name, core[0].name),
                  unpaired=tuple(seg.name for seg in bottom if seg.is_role_token()))
        for side, polarity in (("left", self.left_polarity), ("right", self.right_polarity)):
            top_end = (side == "left") == (polarity == "5-prime")
            sc.overhang("top" if top_end else "bottom",
                        side + " fusion" + ("" if top_end else " bottom"),
                        polarity=polarity, owner="substrate")
        return sc

    def rows(self):
        return self.scene().rows()


@dataclass(frozen=True)
class FivePrimeFragment(CohesiveFragment):
    """Compatibility constructor whose name explicitly fixes both ends as 5′."""

    left_polarity: str = field(default="5-prime", init=False)
    right_polarity: str = field(default="5-prime", init=False)


@dataclass(frozen=True)
class CircularCohesiveAssembly:
    """One source-declared circular ordering of separately supplied fragments.

    IDs are graph input material IDs. The first fragment defines the drawing
    origin; this class never rotates or reverses inputs to find a compatible path.
    """

    fragments: tuple[tuple[str, CohesiveFragment], ...]
    statement: str

    availability = "known"
    action_label = "Present the declared circular fragment order"

    @property
    def ready_label(self):
        polarities = {p for _, f in self.fragments for p in (f.left_polarity, f.right_polarity)}
        return (next(iter(polarities)) if len(polarities) == 1 else "5-prime and 3-prime") + " cohesive fragment ends"

    def __post_init__(self) -> None:
        if len(self.fragments) < 2 or not self.statement.strip():
            raise ValueError("circular assembly needs at least two fragments and a source statement")
        ids = [identifier for identifier, _ in self.fragments]
        if any(not identifier for identifier in ids) or len(set(ids)) != len(ids):
            raise ValueError("assembly material IDs must be nonempty and unique")
        if any(not isinstance(fragment, CohesiveFragment) for _, fragment in self.fragments):
            raise ValueError("every circular assembly input needs explicit cohesive geometry")
        for (_, before), (_, after) in zip(self.fragments, self.fragments[1:] + self.fragments[:1]):
            if (before.right, before.right_polarity) != (after.left, after.left_polarity):
                raise ValueError(f"incompatible cohesive junction: {before.name} → {after.name}")

    def product(self, name: str) -> Construct:
        """Derive one reference strand, recording each shared fusion exactly once."""
        parts = []
        for identifier, fragment in self.fragments:
            parts.append(Segment(identifier + " fusion", fragment.left))
            parts.extend(replace(seg, name=identifier + ": " + seg.name) for seg in fragment.core)
        return Construct(parts, name)

    @property
    def scenes(self) -> tuple[Scene, ...]:
        """Draw nicks and the antiparallel overlap at each source-declared seam."""
        result = []
        for (_, before), (_, after) in zip(self.fragments, self.fragments[1:] + self.fragments[:1]):
            # Only the actual adjacent core segments are needed to show a seam.
            # Role placeholders retain their uncertainty and have no asserted bp length.
            a = replace(before.core[-1], name="upstream core")
            b = replace(after.core[0], name="downstream core")
            ac = complement_segments((a,))[0]
            bc = complement_segments((b,))[0]
            sc = Scene()
            if before.right_polarity == "5-prime":
                sc.strand("upstream top", [a], label=before.name)
                sc.anneal("upstream bottom",
                          [Segment("upstream exposed end", revcomp(before.right)), ac],
                          to="upstream top", pair=(ac.name, a.name),
                          unpaired=(ac.name,) if ac.is_role_token() else ())
                sc.anneal("downstream top", [Segment("downstream exposed end", after.left), b],
                          to="upstream bottom", pair=("downstream exposed end", "upstream exposed end"),
                          label=after.name)
                sc.anneal("downstream bottom", [bc], to="downstream top", pair=(bc.name, b.name),
                          unpaired=(bc.name,) if bc.is_role_token() else ())
            else:
                sc.strand("upstream top", [a, Segment("upstream exposed end", before.right)], label=before.name)
                sc.anneal("upstream bottom", [ac], to="upstream top", pair=(ac.name, a.name),
                          unpaired=(ac.name,) if ac.is_role_token() else ())
                sc.anneal("downstream bottom", [bc, Segment("downstream exposed end", revcomp(after.left))],
                          to="upstream top", pair=("downstream exposed end", "upstream exposed end"), above=True)
                sc.anneal("downstream top", [b], to="downstream bottom", pair=(b.name, bc.name),
                          label=after.name, above=True, unpaired=(b.name,) if b.is_role_token() else ())
            sc.footer("Two nicks remain before sealing; the paired fusion is retained once.",
                      "downstream top" if before.right_polarity == "5-prime" else "downstream bottom",
                      "downstream exposed end")
            result.append(sc)
        return tuple(result)

    def rows(self):
        result = []
        for sc in self.scenes:
            if result:
                result.append(Row())
            result.extend(sc.rows())
        return result
