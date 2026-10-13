"""End repair, dA-tailing and complementary 3'-T adapter ligation."""
from __future__ import annotations

from dataclasses import dataclass

from chemdraw import (ComponentVisual, Construct, Scene, Segment, Workflow, bases_pair,
                      complement_segments, revcomp)
from adapter_ligation import AdapterLigationGeometry


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
                    s.feature, s.length_bp)
            for s in product.insert]
    top_a = Segment("top 3-prime dA", product.three_prime_overhang,
                    inferred=product.overhang_inferred)
    bottom_a = Segment("bottom 3-prime dA", product.three_prime_overhang,
                       inferred=product.overhang_inferred)
    bottom_body = complement_segments(body)
    role_tokens = tuple(s.name for s in bottom_body if s.is_role_token())
    sc = Scene(); sc.strand("top", [*body, top_a], label=label,
                            mod5="p" if product.five_prime_phosphate else "")
    sc.anneal("bottom", [*bottom_body, bottom_a], to="top",
              pair=(bottom_body[-1].name, body[0].name), label=label,
              unpaired=role_tokens, mod5="p" if product.five_prime_phosphate else "")
    sc.overhang("top", top_a.name, polarity="3-prime", owner="substrate", text="3′ dA")
    sc.overhang("bottom", bottom_a.name, polarity="3-prime", owner="substrate", text="3′ dA")
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

    def scene(self) -> Scene:
        """Draw the source-known duplex core and its terminal 3′ overhang."""
        top = [Segment(f"{self.name} duplex core", self.core),
               Segment(f"{self.name} 3-prime overhang", self.overhang_3)]
        bottom = [Segment(f"{self.name} bottom", self.bottom)]
        sc = Scene(); sc.strand("adapter top", top, label=self.name)
        sc.anneal("adapter bottom", bottom, to="adapter top",
                  pair=(bottom[0].name, top[0].name), label=self.name)
        sc.overhang("adapter top", top[-1].name, polarity="3-prime", owner="adapter")
        return sc


def cohesive_ligation_geometry(product: EndPreparedDuplex,
                                adapter: ThreePrimeOverhangAdapter,
                                *, statement: str) -> AdapterLigationGeometry:
    """The two explicit, compatible reagent ends before dA/dT adapter ligation."""
    if not adapter.accepts(product):
        raise ValueError("adapter and end-prepared substrate overhangs are incompatible")
    return AdapterLigationGeometry.known(
        "cohesive", dA_tailed_scene(product), adapter.scene(), statement=statement)


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


def ligation_product_scene(library: Construct, insert_first: str,
                           insert_last: str | None = None, *, label: str = "",
                           adapter_group: str = "adapter",
                           right_adapter_group: str | None = None,
                           highlight_components: bool = True) -> Scene:
    """Draw a two-ended ligation product with component identity made visible.

    Both flanking spans receive the same ``adapter_group`` by construction.  This avoids
    the misleading impression that a symmetric adapter used at both ends represents two
    unrelated reagents.  The insert and both exact covalent boundaries are derived from
    named segments; a missing flank or a non-contiguous boundary cannot be rendered.
    """
    names = [s.name for s in library]
    insert_last = insert_last or insert_first
    try:
        left_i, right_i = names.index(insert_first), names.index(insert_last)
    except ValueError as exc:
        raise ValueError("ligation product insert boundary is not in the library") from exc
    if right_i < left_i:
        raise ValueError("ligation product insert boundaries are reversed")
    if left_i == 0 or right_i == len(names) - 1:
        raise ValueError("two-ended ligation product needs an adapter on both sides")

    right_adapter_group = right_adapter_group or adapter_group
    role_tokens = tuple(s.name + "'" for s in library if s.is_role_token())
    sc = Scene.duplex(list(library), label=label or library.name,
                      unpaired=role_tokens)
    left_last, right_first = names[left_i - 1], names[right_i + 1]
    if highlight_components:
        for strand, suffix in (("top", ""), ("bottom", "'")):
            sc.component(strand, names[0] + suffix, left_last + suffix, adapter_group)
            sc.component(strand, insert_first + suffix, insert_last + suffix, "insert")
            sc.component(strand, right_first + suffix, names[-1] + suffix,
                         right_adapter_group)
    sc.junction("top", left_last, insert_first, "ligation")
    sc.junction("top", insert_last, right_first, "ligation")
    return sc


def react_dA_tail_then_ligate(workflow: Workflow, product: EndPreparedDuplex,
                              adapter: ThreePrimeOverhangAdapter, ligated: Scene, *,
                              dA_action: str = "Add 3-prime dA tails",
                              ligation_action: str = "Ligate adapters",
                              dA_name: str = "dA-tailed DNA",
                              ligated_name: str = "adapter-ligated library",
                              dA_note: str = "",
                              ligation_note: str = "") -> None:
    """Append the two chemically distinct states of sticky-end adapter ligation.

    Callers cannot collapse dA-tailing and ligation into one invisible transition when
    they use this constructor: the ligation-ready overhang state is always rendered.
    """
    ligated_rows = ligated.rows()
    if not any(any(isinstance(v, ComponentVisual) for v in (
            row.visual if isinstance(row.visual, tuple) else (row.visual,)))
            for row in ligated_rows):
        raise ValueError("sticky-end ligation product needs semantic component backgrounds")
    workflow.react(dA_action, dA_tailed_scene(product).rows(), name=dA_name, note=dA_note)
    geometry = cohesive_ligation_geometry(
        product, adapter,
        statement=(f"{adapter.name} carries a 3′-{adapter.overhang_3} overhang "
                   f"complementary to the substrate 3′-{product.three_prime_overhang}."))
    workflow.ligate_adapter(ligation_action, geometry, ligated_rows,
                            name=ligated_name, note=ligation_note)
