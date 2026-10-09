"""End repair, dA-tailing and complementary 3'-T adapter ligation."""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import Construct, Scene, Segment, bases_pair, complement_segments, revcomp


@dataclass(frozen=True)
class EndPreparedDuplex:
    insert: Construct
    five_prime_phosphate: bool = True
    three_prime_overhang: str = "A"
    overhang_inferred: bool = False

    def __post_init__(self) -> None:
        if not len(self.insert):
            raise ValueError("end prep needs a non-empty duplex insert")
        if len(self.three_prime_overhang) != 1:
            raise ValueError("this model represents one-nucleotide 3' overhangs")


def repair_and_dA_tail(insert: Construct, *, inferred: bool = False) -> EndPreparedDuplex:
    """The common polished, 5'-phosphorylated, single-dA-tailed product."""
    return EndPreparedDuplex(insert, five_prime_phosphate=True,
                             three_prime_overhang="A", overhang_inferred=inferred)


def dA_tailed_scene(product: EndPreparedDuplex, label: str = "end-prepared DNA") -> Scene:
    """Draw the two opposed 3' dA overhangs from an :class:`EndPreparedDuplex`."""
    body = [Segment(s.name, s.top, s.tag, s.placeholder, s.inferred, s.bottom, s.note,
                    s.feature)
            for s in product.insert]
    top_a = Segment("top 3-prime dA", product.three_prime_overhang,
                    inferred=product.overhang_inferred)
    bottom_a = Segment("bottom 3-prime dA", product.three_prime_overhang,
                       inferred=product.overhang_inferred)
    bottom_body = complement_segments(body)
    sc = Scene(); sc.strand("top", [*body, top_a], label=label,
                            mod5="p" if product.five_prime_phosphate else "")
    sc.anneal("bottom", [*bottom_body, bottom_a], to="top",
              pair=(bottom_body[0].name, body[0].name), label=label,
              mod5="p" if product.five_prime_phosphate else "")
    sc.mark("top", top_a.name, "3′ dA")
    sc.mark("bottom", bottom_a.name, "3′ dA")
    return sc


@dataclass(frozen=True)
class ThreePrimeOverhangAdapter:
    """A duplex adapter with one extra base at the top strand's 3' end."""
    name: str
    top: str
    bottom: str                 # 5'->3', as synthesized

    def __post_init__(self) -> None:
        if len(self.top) != len(self.bottom) + 1:
            raise ValueError(f"{self.name}: expected a one-base 3' overhang")
        if self.bottom != revcomp(self.top[:-1]):
            raise ValueError(f"{self.name}: duplex core does not base-pair")

    @property
    def core(self) -> str:
        return self.top[:-1]

    @property
    def overhang_3(self) -> str:
        return self.top[-1]

    def accepts(self, product: EndPreparedDuplex) -> bool:
        return (product.five_prime_phosphate
                and bases_pair(self.overhang_3, product.three_prime_overhang))


@dataclass(frozen=True)
class LigationJunction:
    """One sealed adapter/insert boundary, with the paired overhang bases retained."""
    adapter_overhang: str
    insert_overhang: str

    def __post_init__(self) -> None:
        if not bases_pair(self.adapter_overhang, self.insert_overhang):
            raise ValueError("ligation-junction overhangs do not base-pair")

    @property
    def base_pair(self) -> str:
        return f"{self.adapter_overhang}:{self.insert_overhang}"


@dataclass(frozen=True)
class DoubleEndedLigation:
    """A duplex insert sealed to one compatible adapter at each end."""
    insert: EndPreparedDuplex
    adapter: ThreePrimeOverhangAdapter
    left: LigationJunction
    right: LigationJunction

    @property
    def junctions(self) -> tuple[LigationJunction, LigationJunction]:
        return self.left, self.right


def ligate_both_ends(product: EndPreparedDuplex,
                     adapter: ThreePrimeOverhangAdapter) -> DoubleEndedLigation:
    """Construct the sealed two-adapter product, refusing incompatible overhangs."""
    if not adapter.accepts(product):
        raise ValueError(
            f"{adapter.name} {adapter.overhang_3!r} overhang cannot ligate to "
            f"{product.three_prime_overhang!r}-tailed DNA"
        )
    junction = LigationJunction(adapter.overhang_3, product.three_prime_overhang)
    return DoubleEndedLigation(product, adapter, junction, junction)
