"""Reusable structural model for adjacent-probe branched-DNA RNA detection assays."""
from __future__ import annotations

from dataclasses import dataclass
import html

from chemdraw import Scene, Segment, revcomp


@dataclass(frozen=True)
class BranchModel:
    """Published multiplicities for one target RNA and one amplification tree.

    A source may publish the tree only as a total number of label sites, or as
    amplifiers per preamplifier and labels per amplifier.  Give either
    ``labels_per_tree`` alone or both per-level factors (optionally with a matching
    total); a total is never factorized into invented per-level counts.
    """
    target_pairs_min: int
    target_pairs_max: int
    labels_per_tree: int | None = None
    amplifiers_per_preamp: int | None = None
    labels_per_amplifier: int | None = None

    def __post_init__(self) -> None:
        factors = (self.amplifiers_per_preamp, self.labels_per_amplifier)
        if (factors[0] is None) != (factors[1] is None):
            raise ValueError("give both per-level tree factors or neither")
        if factors[0] is not None:
            product = factors[0] * factors[1]
            if self.labels_per_tree is None:
                object.__setattr__(self, "labels_per_tree", product)
            elif self.labels_per_tree != product:
                raise ValueError(f"labels_per_tree {self.labels_per_tree} != "
                                 f"{factors[0]} x {factors[1]}")
        if self.labels_per_tree is None:
            raise ValueError("branched-DNA tree needs labels_per_tree or per-level factors")
        values = (self.target_pairs_min, self.target_pairs_max, self.labels_per_tree,
                  *(f for f in factors if f is not None))
        if any(n <= 0 for n in values) or self.target_pairs_min > self.target_pairs_max:
            raise ValueError(f"invalid branched-DNA multiplicities: {values}")

    @property
    def has_level_factors(self) -> bool:
        return self.amplifiers_per_preamp is not None

    @property
    def labels_per_target(self) -> tuple[int, int]:
        return (self.target_pairs_min * self.labels_per_tree,
                self.target_pairs_max * self.labels_per_tree)


@dataclass(frozen=True)
class SequenceExample:
    """A disclosed implementation without implying identity to a current product."""
    title: str
    note: str
    oligos: tuple[tuple[str, tuple[Segment, ...]], ...]
    scene: Scene
    caption: str


def paired_tail_docking_scene(left_tail: str, right_tail: str,
                              preamplifier_site: str) -> Scene:
    """Draw and validate a preamplifier across two adjacent probe-tail halves.

    The two top segments are separate probe molecules conceptually; joining them in one
    Scene strand supplies a common coordinate system for the 28-base docking duplex.
    The footer makes that representational boundary explicit.
    """
    tails = left_tail + right_tail
    expected = revcomp(tails)
    if preamplifier_site.upper() != expected:
        raise ValueError("preamplifier site is not the reverse complement of both Z tails")
    sc = Scene()
    sc.strand("paired Z tails", [
        Segment("left-probe Z tail", left_tail, "cbc"),
        Segment("right-probe Z tail", right_tail, "cbc"),
    ], label="two adjacent Z-probe tails")
    sc.anneal("preamplifier", [
        Segment("preamplifier docking site", preamplifier_site, "r1"),
    ], to="paired Z tails",
              pair=("preamplifier docking site", "left-probe Z tail"),
              label="preamplifier", above=True)
    sc.mark("paired Z tails", "left-probe Z tail", "tail supplied by probe 1")
    sc.mark("paired Z tails", "right-probe Z tail", "tail supplied by probe 2")
    sc.footer("the two tail segments belong to separate target probes",
              "paired Z tails")
    return sc


def adjacent_target_example_scene(left_target: str, right_target: str) -> Scene:
    """Draw disclosed adjacent target regions with the implied probe-binding arms."""
    sc = Scene()
    sc.strand("target RNA", [
        Segment("left disclosed target region", left_target, "me"),
        Segment("right disclosed target region", right_target, "me"),
    ], label="published target-RNA regions")
    sc.anneal("left probe", [
        Segment("left universal tail", "[undisclosed tail]", "cbc", placeholder=True),
        Segment("left target-binding arm", revcomp(left_target), "r1"),
    ], to="target RNA", pair=("left target-binding arm", "left disclosed target region"),
              label="implied left probe arm", above=True,
              unpaired=("left universal tail",))
    sc.anneal("right probe", [
        Segment("right target-binding arm", revcomp(right_target), "r2"),
        Segment("right universal tail", "[undisclosed tail]", "cbc", placeholder=True),
    ], to="target RNA", pair=("right target-binding arm", "right disclosed target region"),
              label="implied right probe arm", above=True,
              unpaired=("right universal tail",))
    sc.footer("binding arms are derived complements; universal tails were not disclosed",
              "target RNA")
    return sc


def target_scene() -> Scene:
    """One adjacent probe pair, placed by the shared strand-pairing machinery.

    X is an undisclosed base, not a claimed sequence.  The widths are schematic; the
    protocol-specific caption reports any published length range.
    """
    target = [
        Segment("left target site", "X" * 12, placeholder=True),
        Segment("right target site", "X" * 12, placeholder=True),
    ]
    left_probe = [
        Segment("left Z tail", "X" * 6, "cbc", placeholder=True),
        Segment("left probe", "X" * 12, "me", placeholder=True),
    ]
    right_probe = [
        Segment("right probe", "X" * 12, "me", placeholder=True),
        Segment("right Z tail", "X" * 6, "cbc", placeholder=True),
    ]
    sc = Scene()
    sc.strand("target RNA", target, label="target RNA")
    sc.anneal("left Z", left_probe, to="target RNA",
              pair=("left probe", "left target site"),
              label="left target probe", above=True,
              unpaired=("left Z tail",))
    sc.anneal("right Z", right_probe, to="target RNA",
              pair=("right probe", "right target site"),
              label="right target probe", above=True,
              unpaired=("right Z tail",))
    sc.mark("left Z", "left Z tail", "half of preamplifier landing site")
    sc.mark("right Z", "right Z tail", "half of preamplifier landing site")
    sc.footer("schematic widths; X = undisclosed base", "target RNA")
    return sc


def tree_svg(model: BranchModel) -> str:
    """Render the amplification hierarchy as a molecular graph, not ASCII punctuation."""
    width, height = 860, 430
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" class="chem-svg topology-svg" '
             f'width="{width}" height="{height}" viewBox="0 0 {width} {height}" '
             'role="img" aria-label="branched-DNA signal-amplification topology">']
    parts += [
        '<line x1="90" y1="386" x2="770" y2="386" stroke="var(--ink)" stroke-width="4"/>',
        '<text x="82" y="391" text-anchor="end" class="topology-end">5′</text>',
        '<text x="778" y="391" class="topology-end">3′ target RNA</text>',
        '<path d="M350,386 L350,350 L420,322" fill="none" stroke="var(--c-r1)" stroke-width="6"/>',
        '<path d="M510,386 L510,350 L440,322" fill="none" stroke="var(--c-r2)" stroke-width="6"/>',
        '<text x="430" y="410" text-anchor="middle" class="topology-comment">one adjacent target-probe pair</text>',
        '<line x1="430" y1="322" x2="430" y2="235" stroke="var(--c-tso)" stroke-width="8"/>',
        '<text x="445" y="277" class="topology-label">preamplifier trunk</text>',
    ]
    # Five representative amplifier branches show topology; the label states the exact
    # model multiplicity instead of drawing dozens of indistinguishable arms.
    amplifier_x = (190, 310, 430, 550, 670)
    for x in amplifier_x:
        parts.append(f'<path d="M430,235 Q430,205 {x},188 L{x},132" fill="none" '
                     'stroke="var(--c-cbc)" stroke-width="6"/>')
        for dx in (-28, -14, 0, 14, 28):
            parts.append(f'<line x1="{x}" y1="145" x2="{x+dx}" y2="88" '
                         'stroke="var(--c-umi)" stroke-width="3"/>')
            parts.append(f'<circle cx="{x+dx}" cy="82" r="5" fill="var(--c-umi)"/>')
    if model.has_level_factors:
        parts += [
            f'<text x="430" y="215" text-anchor="middle" class="topology-label">'
            f'{model.amplifiers_per_preamp} amplifier branches per preamplifier</text>',
            f'<text x="430" y="54" text-anchor="middle" class="topology-label">'
            f'{model.labels_per_amplifier} label probes per amplifier</text>',
        ]
    else:
        parts += [
            '<text x="430" y="215" text-anchor="middle" class="topology-label">'
            'amplifier branches (per-level count not modelled)</text>',
            '<text x="430" y="54" text-anchor="middle" class="topology-label">'
            'label probes (per-amplifier count not modelled)</text>',
        ]
    parts += [
        f'<text x="430" y="26" text-anchor="middle" class="topology-heading">'
        f'{model.labels_per_tree} potential labels per signal tree</text>',
        '</svg>',
    ]
    return "".join(parts)


def tree_panel(model: BranchModel, caption: str = "") -> str:
    cap = f"<figcaption>{html.escape(caption)}</figcaption>" if caption else ""
    return (f'<figure class="chem-panel topology-panel">{cap}'
            f'<div class="diagram-scroll">{tree_svg(model)}</div></figure>')
