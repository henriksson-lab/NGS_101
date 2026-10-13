"""Constraint-checked topological drawings for non-linear nucleic-acid assemblies.

Sequence rows are ideal for linear libraries, but they obscure loops, circles and
multi-strand scaffolds.  This module keeps molecular identity and pairing as data while
using a schematic path layout.  It does not predict thermodynamic folding: callers state
which spans are intended to pair, and constructors reject incompatible sequences.

The first reusable topology is a two-probe proximity-templated circle (PLAYR and related
assays).  The underlying SVG vocabulary -- semantic strand paths, paired spans and sealed
junctions -- is deliberately independent of RNA versus DNA.
"""
from __future__ import annotations

from dataclasses import dataclass
import html
import math

from chemdraw import Segment, revcomp


Point = tuple[float, float]


def _sample_polyline(vertices: list[Point], count: int) -> list[Point]:
    """Place ``count`` nucleotide centres uniformly along a routed polyline."""
    if count <= 0:
        return []
    lengths = [math.dist(a, b) for a, b in zip(vertices, vertices[1:])]
    total = sum(lengths)
    if total == 0:
        raise ValueError("a flexible route must have non-zero length")
    out = []
    for i in range(1, count + 1):
        wanted = total * i / (count + 1)
        walked = 0.0
        for (a, b), length in zip(zip(vertices, vertices[1:]), lengths):
            if wanted <= walked + length:
                fraction = (wanted - walked) / length
                out.append((a[0] + (b[0] - a[0]) * fraction,
                            a[1] + (b[1] - a[1]) * fraction))
                break
            walked += length
    return out


def _orientation(a: Point, b: Point, c: Point) -> float:
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])


def _proper_cross(a: Point, b: Point, c: Point, d: Point) -> bool:
    """True for a crossing, overlap or unintended touch; shared endpoints are joints."""
    if a in (c, d) or b in (c, d):
        return False
    ab_c, ab_d = _orientation(a, b, c), _orientation(a, b, d)
    cd_a, cd_b = _orientation(c, d, a), _orientation(c, d, b)
    if ab_c * ab_d < -1e-6 and cd_a * cd_b < -1e-6:
        return True

    def on_segment(p: Point, q: Point, r: Point) -> bool:
        return (abs(_orientation(p, q, r)) <= 1e-6
                and min(p[0], q[0]) - 1e-6 <= r[0] <= max(p[0], q[0]) + 1e-6
                and min(p[1], q[1]) - 1e-6 <= r[1] <= max(p[1], q[1]) + 1e-6)
    return any((on_segment(a, b, c), on_segment(a, b, d),
                on_segment(c, d, a), on_segment(c, d, b)))


def _route_cost(vertices: list[Point]) -> float:
    """Aesthetic energy: compact paths with few, gentle changes of direction."""
    cost = sum(math.dist(a, b) for a, b in zip(vertices, vertices[1:]))
    for a, b, c in zip(vertices, vertices[1:], vertices[2:]):
        incoming = (b[0] - a[0], b[1] - a[1])
        outgoing = (c[0] - b[0], c[1] - b[1])
        ilen = math.hypot(*incoming); olen = math.hypot(*outgoing)
        cosine = (incoming[0] * outgoing[0] + incoming[1] * outgoing[1]) / (ilen * olen)
        cost += 70 + 35 * (1 - cosine)
    return cost


def _choose_route(candidates: list[list[Point]], blockers: list[tuple[Point, Point]]) -> list[Point]:
    """Choose the least-complex non-crossing route from deterministic candidates."""
    valid = []
    for route in candidates:
        edges = list(zip(route, route[1:]))
        if any(_proper_cross(a, b, c, d) for a, b in edges for c, d in blockers):
            continue
        if any(_proper_cross(a, b, c, d) for i, (a, b) in enumerate(edges)
               for c, d in edges[i + 2:]):
            continue
        valid.append(route)
    if not valid:
        raise ValueError("no crossing-free route satisfies the fold hints")
    return min(valid, key=_route_cost)


@dataclass(frozen=True)
class FoldStrand:
    """A nucleotide-resolved strand after layout, in molecular 5′→3′ order."""

    name: str
    sequence: str
    points: tuple[Point, ...]
    colour: str = "var(--ink)"
    polymer: str = "DNA"
    label: str = ""

    def __post_init__(self) -> None:
        sequence = _clean(self.sequence, self.name)
        object.__setattr__(self, "sequence", sequence)
        if len(sequence) != len(self.points):
            raise ValueError(f"{self.name}: one point is required for every nucleotide")


@dataclass(frozen=True)
class FoldPairing:
    """An antiparallel, base-resolved pairing constraint between two strands."""

    first: str
    first_start: int
    second: str
    second_start: int
    length: int
    name: str


@dataclass(frozen=True)
class StraightHint:
    """Require a strand interval to remain collinear in the final embedding."""

    strand: str
    start: int
    end: int


class FoldLayout:
    """Validated planar fold; layout is geometric, never a free-energy prediction."""

    def __init__(self, strands: list[FoldStrand], pairings: list[FoldPairing],
                 extra_bonds: list[tuple[str, int, str, int, str]] = (),
                 straight: list[StraightHint] = ()) -> None:
        self.strands = {strand.name: strand for strand in strands}
        if len(self.strands) != len(strands):
            raise ValueError("fold strand names must be unique")
        self.pairings = tuple(pairings)
        self.extra_bonds = tuple(extra_bonds)
        self.straight = tuple(straight)
        self._validate_pairing()
        self._validate_hints()
        self._validate_no_crossings()

    def flexible_runs(self) -> tuple[tuple[str, int, int], ...]:
        """Infer maximal unpaired intervals; these are the layout's flexible joints."""
        paired = {name: set() for name in self.strands}
        for pairing in self.pairings:
            for a, ai, b, bi in self._paired_indices(pairing):
                paired[a.name].add(ai); paired[b.name].add(bi)
        runs = []
        for name, strand in self.strands.items():
            start = None
            for index in range(len(strand.sequence) + 1):
                is_free = index < len(strand.sequence) and index not in paired[name]
                if is_free and start is None:
                    start = index
                elif not is_free and start is not None:
                    runs.append((name, start, index - 1)); start = None
        return tuple(runs)

    def _paired_indices(self, pairing: FoldPairing):
        a, b = self.strands[pairing.first], self.strands[pairing.second]
        for offset in range(pairing.length):
            yield a, pairing.first_start + offset, b, pairing.second_start + pairing.length - 1 - offset

    def _validate_pairing(self) -> None:
        for pairing in self.pairings:
            if pairing.first not in self.strands or pairing.second not in self.strands:
                raise ValueError(f"{pairing.name}: pairing names an unknown strand")
            for a, ai, b, bi in self._paired_indices(pairing):
                if not (0 <= ai < len(a.sequence) and 0 <= bi < len(b.sequence)):
                    raise ValueError(f"{pairing.name}: pairing extends beyond a strand")
                if a.sequence[ai] != "N" and b.sequence[bi] != "N" \
                        and a.sequence[ai] != revcomp(b.sequence[bi]):
                    raise ValueError(f"{pairing.name}: {a.sequence[ai]} does not pair with {b.sequence[bi]}")

    def _validate_hints(self) -> None:
        for hint in self.straight:
            if hint.strand not in self.strands:
                raise ValueError("straight hint names an unknown strand")
            points = self.strands[hint.strand].points[hint.start:hint.end + 1]
            if len(points) < 2 or hint.start < 0 or hint.end >= len(self.strands[hint.strand].points):
                raise ValueError("straight hint has an invalid interval")
            if any(abs(_orientation(points[0], points[-1], point)) > 1e-6
                   for point in points[1:-1]):
                raise ValueError(f"straight hint was not satisfied for {hint.strand}")

    def _edges(self):
        edges = []
        for strand in self.strands.values():
            edges.extend((a, b, f"{strand.name} backbone")
                         for a, b in zip(strand.points, strand.points[1:]))
        for pairing in self.pairings:
            edges.extend((a.points[ai], b.points[bi], f"{pairing.name} base pair")
                         for a, ai, b, bi in self._paired_indices(pairing))
        for first, fi, second, si, name in self.extra_bonds:
            edges.append((self.strands[first].points[fi], self.strands[second].points[si], name))
        return edges

    def _validate_no_crossings(self) -> None:
        edges = self._edges()
        for i, (a, b, aname) in enumerate(edges):
            for c, d, bname in edges[i + 1:]:
                if _proper_cross(a, b, c, d):
                    raise ValueError(f"fold is not planar: {aname} crosses {bname}")

    @staticmethod
    def _path(points: tuple[Point, ...], colour: str) -> str:
        coords = " ".join(f"{x:.1f},{y:.1f}" for x, y in points)
        return (f'<polyline points="{coords}" fill="none" stroke="{colour}" '
                'stroke-width="2" stroke-linejoin="round"/>')

    def svg(self, width: int, height: int, title: str) -> str:
        parts = [f'<svg xmlns="http://www.w3.org/2000/svg" class="chem-svg fold-svg" '
                 f'width="{width}" height="{height}" viewBox="0 0 {width} {height}" '
                 f'role="img" aria-label="{html.escape(title, quote=True)}">']
        # Pairing rungs sit behind the strands and letters.
        for pairing in self.pairings:
            for a, ai, b, bi in self._paired_indices(pairing):
                x1, y1 = a.points[ai]; x2, y2 = b.points[bi]
                parts.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
                             'class="fold-pair"><title>' + html.escape(pairing.name) + '</title></line>')
        closed_strands = {
            first for first, fi, second, si, _ in self.extra_bonds
            if first == second and {fi, si} == {0, len(self.strands[first].sequence) - 1}
        }
        for strand in self.strands.values():
            parts.append(f'<g class="fold-strand"><title>{html.escape(Segment(strand.name, strand.sequence).hover_text())}</title>')
            parts.append(self._path(strand.points, strand.colour))
            for base, (x, y) in zip(strand.sequence, strand.points):
                parts.append(f'<text x="{x:.1f}" y="{y + 4:.1f}" text-anchor="middle" '
                             f'fill="{strand.colour}" class="fold-base">{base}</text>')
            if strand.name not in closed_strands:
                x5, y5 = strand.points[0]; x3, y3 = strand.points[-1]
                parts.append(f'<text x="{x5:.1f}" y="{y5 - 11:.1f}" text-anchor="middle" class="fold-end">5′</text>')
                parts.append(f'<text x="{x3:.1f}" y="{y3 - 11:.1f}" text-anchor="middle" class="fold-end">3′</text>')
            parts.append('</g>')
        for first, fi, second, si, name in self.extra_bonds:
            x1, y1 = self.strands[first].points[fi]; x2, y2 = self.strands[second].points[si]
            parts.append(f'<g><title>{html.escape(name)}</title><line x1="{x1:.1f}" y1="{y1:.1f}" '
                         f'x2="{x2:.1f}" y2="{y2:.1f}" class="fold-ligation"/>'
                         f'<circle cx="{(x1+x2)/2:.1f}" cy="{(y1+y2)/2:.1f}" r="4" class="topology-junction"/></g>')
        legend_x = 18.0
        for strand in self.strands.values():
            label = strand.label or strand.name
            parts.append(f'<line x1="{legend_x:.1f}" y1="{height-18}" x2="{legend_x+16:.1f}" '
                         f'y2="{height-18}" stroke="{strand.colour}" stroke-width="4"/>'
                         f'<text x="{legend_x+22:.1f}" y="{height-14}" class="fold-label">'
                         f'{html.escape(label)}</text>')
            legend_x += 34 + len(label) * 7.2
        parts.append('</svg>')
        return "".join(parts)


def _clean(sequence: str, name: str) -> str:
    value = sequence.upper().replace(" ", "")
    if not value or set(value) - set("ACGTUN"):
        raise ValueError(f"{name} must be a non-empty nucleic-acid sequence")
    return value


@dataclass(frozen=True)
class ProximityProbe:
    """One 5′ target arm followed by a linker and two 3′ docking spans."""

    name: str
    target_arm: str
    linker: str
    first_role: str
    first_site: str
    second_role: str
    second_site: str

    def __post_init__(self) -> None:
        for attr in ("target_arm", "linker", "first_site", "second_site"):
            object.__setattr__(self, attr, _clean(getattr(self, attr), f"{self.name} {attr}"))
        if {self.first_role, self.second_role} != {"inner", "outer"}:
            raise ValueError("a proximity probe must expose one inner and one outer site")

    def site(self, role: str) -> str:
        if self.first_role == role:
            return self.first_site
        if self.second_role == role:
            return self.second_site
        raise KeyError(role)

    @property
    def sequence(self) -> str:
        return self.target_arm + self.linker + self.first_site + self.second_site


@dataclass(frozen=True)
class ProximityCircleTopology:
    """Two adjacent target probes that template an inner bridge and outer loop.

    ``inner`` is written 5′→3′ from the downstream probe to the upstream probe.
    ``outer`` is written 5′→3′ from the upstream probe around the loop to the downstream
    probe.  Consequently their ends meet at exactly two ligation junctions.
    """

    upstream: ProximityProbe
    downstream: ProximityProbe
    inner: str
    outer: str
    inner_spacer: str
    outer_body: str
    name: str = "proximity-templated reporter circle"

    def __post_init__(self) -> None:
        for attr in ("inner", "outer", "inner_spacer", "outer_body"):
            object.__setattr__(self, attr, _clean(getattr(self, attr), attr))
        if (self.upstream.first_role, self.upstream.second_role) != ("outer", "inner"):
            raise ValueError("upstream probe sites must run outer→inner toward the junction")
        if (self.downstream.first_role, self.downstream.second_role) != ("inner", "outer"):
            raise ValueError("downstream probe sites must run inner→outer away from the junction")
        expected_inner = (revcomp(self.downstream.site("inner")) + self.inner_spacer
                          + revcomp(self.upstream.site("inner")))
        expected_outer = (revcomp(self.upstream.site("outer")) + self.outer_body
                          + revcomp(self.downstream.site("outer")))
        if self.inner != expected_inner:
            raise ValueError("inner oligo does not bridge the declared probe sites")
        if self.outer != expected_outer:
            raise ValueError("outer oligo does not bridge the declared probe sites")

    @staticmethod
    def _title(name: str, sequence: str) -> str:
        return html.escape(Segment(name, sequence).hover_text(), quote=True)

    @staticmethod
    def _line(x1, y1, x2, y2, colour, title, width=7, extra="") -> str:
        return (f'<g><title>{title}</title><path d="M{x1},{y1} L{x2},{y2}" '
                f'stroke="{colour}" stroke-width="{width}" stroke-linecap="round" '
                f'fill="none" {extra}/></g>')

    def _probe(self, x: int, y_rna: int, probe: ProximityProbe) -> str:
        scale = 5
        target_w = max(92, len(probe.target_arm) * 3)
        site_w = len(probe.first_site) * scale
        site2_w = len(probe.second_site) * scale
        y_top = y_rna - 92
        x0, x1 = x, x + target_w
        # Molecular path: 5′ at target-arm right, antiparallel along RNA, then the
        # 3′ linker turns upward and the two docking sites run left-to-right.
        target_title = self._title(f"{probe.name} RNA-binding arm", probe.target_arm)
        linker_title = self._title(f"{probe.name} poly(A) linker", probe.linker)
        first_title = self._title(f"{probe.name} {probe.first_role} docking site",
                                  probe.first_site)
        second_title = self._title(f"{probe.name} {probe.second_role} docking site",
                                   probe.second_site)
        colour = {"inner": "var(--c-w1)", "outer": "var(--c-s5)"}
        parts = [
            self._line(x1, y_rna - 10, x0, y_rna - 10, "var(--c-r1)", target_title),
            self._line(x0, y_rna - 10, x0, y_top, "var(--c-tso)", linker_title),
            self._line(x0, y_top, x0 + site_w, y_top, colour[probe.first_role], first_title),
            self._line(x0 + site_w, y_top, x0 + site_w + site2_w, y_top,
                       colour[probe.second_role], second_title),
            f'<text x="{x0 + target_w / 2}" y="{y_rna - 20}" text-anchor="middle" '
            f'class="topology-label">{html.escape(probe.name)}</text>',
            f'<text x="{x1 + 5}" y="{y_rna - 5}" class="topology-end">5′</text>',
            f'<text x="{x0 + site_w + site2_w + 5}" y="{y_top + 5}" '
            f'class="topology-end">3′</text>',
        ]
        return "".join(parts)

    def _stage(self, xoff: int, with_circle: bool) -> str:
        y_rna = 245
        up_x, down_x = xoff + 55, xoff + 260
        parts = [
            self._line(xoff + 20, y_rna, xoff + 430, y_rna, "var(--ink)",
                       html.escape("target RNA; adjacent probe-binding sites"), width=4),
            f'<text x="{xoff + 10}" y="{y_rna + 6}" class="topology-end">5′</text>',
            f'<text x="{xoff + 437}" y="{y_rna + 6}" class="topology-end">3′ RNA</text>',
            self._probe(up_x, y_rna, self.upstream),
            self._probe(down_x, y_rna, self.downstream),
        ]
        if not with_circle:
            parts.append(f'<text x="{xoff + 225}" y="300" text-anchor="middle" '
                         'class="topology-comment">adjacent 5′ target arms; exposed 3′ docking platform</text>')
            return "".join(parts)

        scale = 5
        up_first = len(self.upstream.first_site) * scale
        # Upstream is outer then inner; downstream is inner then outer for PLAYR-like
        # assemblies.  The constructor permits either order, so calculate the spans.
        def spans(x, probe):
            a = len(probe.first_site) * scale
            b = len(probe.second_site) * scale
            return ({probe.first_role: (x, x + a),
                     probe.second_role: (x + a, x + a + b)})
        us, ds = spans(up_x, self.upstream), spans(down_x, self.downstream)
        y_pair = y_rna - 102
        # Inner oligo: 5′ on downstream, 3′ on upstream, with its unpaired spacer
        # routed through the gap.  The visual path follows molecular direction.
        inner_title = self._title("phosphorylated inner bridge", self.inner)
        ix1, ix2 = ds["inner"][1], us["inner"][0]
        mid = (ix1 + ix2) / 2
        parts.append(
            f'<g><title>{inner_title}</title><path d="M{ix1},{y_pair} '
            f'L{ds["inner"][0]},{y_pair} Q{mid},{y_pair + 26} {us["inner"][1]},{y_pair} '
            f'L{us["inner"][0]},{y_pair}" stroke="var(--c-cbc)" stroke-width="6" '
            'stroke-linecap="round" fill="none"/></g>')
        parts.append(f'<text x="{mid}" y="{y_pair + 35}" text-anchor="middle" '
                     'class="topology-label">inner / insert</text>')
        parts.append(f'<text x="{ix1 + 5}" y="{y_pair - 7}" class="topology-end">5′</text>')
        parts.append(f'<text x="{us["inner"][0] - 18}" y="{y_pair - 7}" class="topology-end">3′</text>')

        # Outer oligo: 5′ on upstream, around its long body, 3′ on downstream.
        outer_title = self._title("phosphorylated outer loop", self.outer)
        # Both docking platforms run 5′→3′ left-to-right, so the complementary
        # outer strand traverses each paired span right-to-left.
        up_outer_left, up_outer_right = us["outer"]
        down_outer_left, down_outer_right = ds["outer"]
        y_outer = 55
        parts.append(
            f'<g><title>{outer_title}</title><path d="M{up_outer_right},{y_pair} '
            f'L{up_outer_left},{y_pair} Q{up_outer_left - 18},{y_outer} '
            f'{up_outer_left + 55},{y_outer} L{down_outer_right - 55},{y_outer} '
            f'Q{down_outer_right + 18},{y_outer} {down_outer_right},{y_pair} '
            f'L{down_outer_left},{y_pair}" stroke="var(--note-rule)" stroke-width="6" '
            'stroke-linecap="round" fill="none"/></g>')
        parts.append(f'<text x="{(up_outer_left + down_outer_right)/2}" y="{y_outer - 10}" '
                     'text-anchor="middle" class="topology-label">outer / backbone</text>')
        parts.append(f'<text x="{up_outer_right + 5}" y="{y_pair - 7}" class="topology-end">5′</text>')
        parts.append(f'<text x="{down_outer_left - 18}" y="{y_pair - 7}" class="topology-end">3′</text>')
        # The two coloured dots are the only covalent closures.
        for x in (up_outer_right, ds["inner"][1]):
            parts.append(f'<circle cx="{x}" cy="{y_pair}" r="6" class="topology-junction">'
                         '<title>T4 DNA ligase junction</title></circle>')
        parts.append(f'<text x="{xoff + 225}" y="300" text-anchor="middle" '
                     'class="topology-comment">two 5′-phosphorylated oligos; two sealed junctions</text>')
        return "".join(parts)

    def svg(self) -> str:
        """Render probe binding and circle assembly side by side."""
        return (f'<svg xmlns="http://www.w3.org/2000/svg" class="chem-svg topology-svg" '
                f'width="980" height="330" viewBox="0 0 980 330" role="img" '
                f'aria-label="{html.escape(self.name, quote=True)}">'
                '<text x="235" y="25" text-anchor="middle" class="topology-heading">'
                '1. Adjacent probe pair on target</text>'
                '<text x="735" y="25" text-anchor="middle" class="topology-heading">'
                '2. Insert and backbone hybridize</text>'
                + self._stage(0, False)
                + '<path d="M468,170 L500,170" class="topology-process"/>'
                  '<path d="M500,170 l-10,-7 v14 z" class="topology-process-fill"/>'
                + self._stage(500, True) + '</svg>')

    def folded_layout(self) -> FoldLayout:
        """Resolve the complete ligation-ready complex at single-base resolution.

        The sole positional hint is a straight target RNA. Pairing constraints place the
        four helices. Runs not covered by a pairing constraint are flexible joints; the
        outer body is routed through the lowest-complexity crossing-free lane.
        """
        step = 12.0
        target_points = tuple((152 + step * i, 380.0) for i in range(48))
        target = FoldStrand("target RNA (symbolic adjacent sites)", "N" * 48,
                            target_points, "var(--ink)", "RNA", "target RNA")

        def probe_points(rna_points: tuple[Point, ...], site_x: float, probe: ProximityProbe):
            arm = tuple((x, 356.0) for x, _ in reversed(rna_points))
            joint_start = arm[-1]
            linker_route = [joint_start, (site_x - 24, 328), (site_x - 24, 220),
                            (site_x, 180)]
            linker = tuple(_sample_polyline(linker_route, len(probe.linker)))
            sites = tuple((site_x + step * i, 180.0)
                          for i in range(len(probe.first_site) + len(probe.second_site)))
            return arm + linker + sites

        upstream_points = probe_points(target_points[:24], 164, self.upstream)
        downstream_points = probe_points(target_points[24:], 440, self.downstream)
        upstream = FoldStrand(self.upstream.name, self.upstream.sequence, upstream_points,
                              "var(--c-r1)", label="Probe 2")
        downstream = FoldStrand(self.downstream.name, self.downstream.sequence,
                                downstream_points, "var(--c-r2)", label="Probe 1")

        # The paired portions are fixed by the probe sites. The three-base insert spacer
        # is the direct (zero-bend) solution between them.
        insert_first = tuple((536 - step * i, 156.0) for i in range(9))
        insert_last = tuple((368 - step * i, 156.0) for i in range(9))
        insert_joint = tuple(_sample_polyline([insert_first[-1], insert_last[0]],
                                              len(self.inner_spacer)))
        insert = FoldStrand("phosphorylated insert", self.inner,
                            insert_first + insert_joint + insert_last, "var(--c-cbc)",
                            label="insert")

        outer_first = tuple((260 - step * i, 156.0) for i in range(9))
        outer_last = tuple((644 - step * i, 156.0) for i in range(9))
        blockers = [((152, 180), (644, 180)),
                    ((insert_last[-1][0], 156), (insert_first[0][0], 156))]
        routes = [[outer_first[-1], (140, 132), (140, lane), (668, lane),
                   (668, 132), outer_last[0]] for lane in (108, 84, 60)]
        outer_route = _choose_route(routes, blockers)
        outer_joint = tuple(_sample_polyline(outer_route, len(self.outer_body)))
        outer = FoldStrand("phosphorylated universal backbone", self.outer,
                           outer_first + outer_joint + outer_last, "var(--c-s5)",
                           label="backbone")

        probe_site = len(self.upstream.target_arm) + len(self.upstream.linker)
        pairings = [
            FoldPairing(target.name, 0, upstream.name, 0, 24, "Probe 2–target RNA"),
            FoldPairing(target.name, 24, downstream.name, 0, 24, "Probe 1–target RNA"),
            FoldPairing(upstream.name, probe_site, outer.name, 0, 9,
                        "Probe 2 outer/backbone site"),
            FoldPairing(upstream.name, probe_site + 9, insert.name, 12, 9,
                        "Probe 2 inner/insert site"),
            FoldPairing(downstream.name, probe_site, insert.name, 0, 9,
                        "Probe 1 inner/insert site"),
            FoldPairing(downstream.name, probe_site + 9, outer.name, len(self.outer) - 9, 9,
                        "Probe 1 outer/backbone site"),
        ]
        junctions = [(insert.name, len(self.inner) - 1, outer.name, 0,
                      "T4 DNA ligase junction"),
                     (outer.name, len(self.outer) - 1, insert.name, 0,
                      "T4 DNA ligase junction")]
        return FoldLayout([target, upstream, downstream, insert, outer], pairings, junctions,
                          [StraightHint(target.name, 0, len(target.sequence) - 1)])

    def folded_svg(self) -> str:
        return self.folded_layout().svg(820, 450, self.name + ", nucleotide-resolved fold")

    def folded_panel(self, caption: str = "") -> str:
        cap = f"<figcaption>{html.escape(caption)}</figcaption>" if caption else ""
        return (f'<figure class="chem-panel topology-panel">{cap}'
                f'<div class="diagram-scroll">{self.folded_svg()}</div></figure>')
    def panel(self, caption: str = "") -> str:
        cap = f"<figcaption>{html.escape(caption)}</figcaption>" if caption else ""
        return (f'<figure class="chem-panel topology-panel">{cap}'
                f'<div class="diagram-scroll">{self.svg()}</div></figure>')


@dataclass(frozen=True)
class HairpinTopology:
    """A sequence-resolved single strand with two declared antiparallel stem spans."""

    name: str
    sequence: str
    first_stem_start: int
    second_stem_start: int
    stem_length: int
    colour: str = "var(--c-r1)"

    def __post_init__(self) -> None:
        sequence = _clean(self.sequence, self.name)
        object.__setattr__(self, "sequence", sequence)
        a, b, length = self.first_stem_start, self.second_stem_start, self.stem_length
        if length < 2 or not (0 <= a < a + length <= b < b + length <= len(sequence)):
            raise ValueError("hairpin stem spans must be ordered, disjoint and in range")
        if sequence[a:a + length] != revcomp(sequence[b:b + length]):
            raise ValueError("declared hairpin stem spans are not reverse-complementary")

    @property
    def loop_start(self) -> int:
        return self.first_stem_start + self.stem_length

    @property
    def loop_length(self) -> int:
        return self.second_stem_start - self.loop_start

    def folded_layout(self) -> FoldLayout:
        step = 12.0
        prefix = self.first_stem_start
        suffix = len(self.sequence) - (self.second_stem_start + self.stem_length)
        x0 = 60 + max(prefix, suffix) * step
        y_first, y_second = 210.0, 170.0
        first = tuple((x0 + step * i, y_first) for i in range(self.stem_length))
        second = tuple((x0 + step * (self.stem_length - 1 - i), y_second)
                       for i in range(self.stem_length))
        prefix_points = tuple((x0 - step * (prefix - i), y_first) for i in range(prefix))
        suffix_points = tuple((x0 - step * (i + 1), y_second) for i in range(suffix))
        extent = max(32.0, self.loop_length * step / 3)
        right = first[-1][0]
        loop_route = [first[-1], (right + extent, y_first),
                      (right + extent, y_second), second[0]]
        loop = tuple(_sample_polyline(loop_route, self.loop_length))
        points = prefix_points + first + loop + second + suffix_points
        strand = FoldStrand(self.name, self.sequence, points, self.colour,
                            label=self.name)
        pairing = FoldPairing(self.name, self.first_stem_start, self.name,
                              self.second_stem_start, self.stem_length,
                              f"{self.name} stem")
        straight = [StraightHint(self.name, 0, self.loop_start - 1),
                    StraightHint(self.name, self.second_stem_start, len(self.sequence) - 1)]
        return FoldLayout([strand], [pairing], straight=straight)

    def folded_svg(self) -> str:
        prefix = self.first_stem_start
        suffix = len(self.sequence) - (self.second_stem_start + self.stem_length)
        width = int(150 + max(prefix, suffix) * 12 + self.stem_length * 12
                    + max(32, self.loop_length * 4))
        return self.folded_layout().svg(max(520, width), 280, self.name + " hairpin fold")

    def folded_panel(self, caption: str = "") -> str:
        cap = f"<figcaption>{html.escape(caption)}</figcaption>" if caption else ""
        return (f'<figure class="chem-panel topology-panel">{cap}'
                f'<div class="diagram-scroll">{self.folded_svg()}</div></figure>')


def _closed_rectangle_points(count: int, step: float = 12.0) -> tuple[tuple[Point, ...], int, int]:
    """Lay a covalent circle on a compact, orthogonal, non-self-crossing perimeter."""
    if count < 4:
        raise ValueError("a drawn nucleotide circle needs at least four bases")
    width = min(600.0, max(144.0, count * step * 0.38))
    height = max(88.0, (count * step - 2 * width) / 2)
    x0, y0 = 64.0, 58.0
    vertices = [(x0, y0), (x0 + width, y0), (x0 + width, y0 + height),
                (x0, y0 + height), (x0, y0)]
    lengths = [math.dist(a, b) for a, b in zip(vertices, vertices[1:])]
    perimeter = sum(lengths)
    points = []
    for index in range(count):
        wanted = perimeter * index / count
        walked = 0.0
        for (a, b), length in zip(zip(vertices, vertices[1:]), lengths):
            if wanted <= walked + length:
                fraction = (wanted - walked) / length
                points.append((a[0] + (b[0] - a[0]) * fraction,
                               a[1] + (b[1] - a[1]) * fraction))
                break
            walked += length
    return tuple(points), int(x0 + width + 70), int(y0 + height + 92)


@dataclass(frozen=True)
class ClosedCircleTopology:
    """A nucleotide-resolved covalently closed single strand."""

    name: str
    sequence: str
    colour: str = "var(--c-r3)"

    def __post_init__(self) -> None:
        object.__setattr__(self, "sequence", _clean(self.sequence, self.name))

    def folded_layout(self) -> tuple[FoldLayout, int, int]:
        points, width, height = _closed_rectangle_points(len(self.sequence))
        strand = FoldStrand(self.name, self.sequence, points, self.colour,
                            label="closed template")
        closure = [(strand.name, len(strand.sequence) - 1, strand.name, 0,
                    "covalently closed circle")]
        return FoldLayout([strand], [], closure), width, height

    def folded_svg(self) -> str:
        layout, width, height = self.folded_layout()
        return layout.svg(width, height, self.name)

    def folded_panel(self, caption: str = "") -> str:
        cap = f"<figcaption>{html.escape(caption)}</figcaption>" if caption else ""
        return (f'<figure class="chem-panel topology-panel">{cap}'
                f'<div class="diagram-scroll">{self.folded_svg()}</div></figure>')


@dataclass(frozen=True)
class CircularPrimerTopology:
    """A closed ssDNA template with one exact antiparallel primer-binding site."""

    name: str
    circle_sequence: str
    primer_name: str
    primer_sequence: str
    circle_colour: str = "var(--c-r3)"
    primer_colour: str = "var(--c-r2)"

    def __post_init__(self) -> None:
        circle = _clean(self.circle_sequence, self.name)
        primer = _clean(self.primer_sequence, self.primer_name)
        object.__setattr__(self, "circle_sequence", circle)
        object.__setattr__(self, "primer_sequence", primer)
        if len(primer) >= len(circle):
            raise ValueError("a circle-bound primer must be shorter than its template")
        site = revcomp(primer)
        doubled = circle + circle[:len(site) - 1]
        hits = [index for index in range(len(circle)) if doubled.startswith(site, index)]
        if len(hits) != 1:
            raise ValueError(f"{self.primer_name} has {len(hits)} circular binding sites")

    def folded_layout(self) -> tuple[FoldLayout, int, int]:
        site = revcomp(self.primer_sequence)
        doubled = self.circle_sequence + self.circle_sequence[:len(site) - 1]
        start = next(index for index in range(len(self.circle_sequence))
                     if doubled.startswith(site, index))
        rotated = self.circle_sequence[start:] + self.circle_sequence[:start]
        perimeter_points, width, height = _closed_rectangle_points(len(rotated))
        top_count = next((i for i, (_, y) in enumerate(perimeter_points)
                          if abs(y - perimeter_points[0][1]) > 1e-6), len(perimeter_points))
        spare = top_count - len(self.primer_sequence)
        if spare < 2:
            raise ValueError("circular primer site needs space on both sides of one edge")
        margin = spare // 2
        # Sequence position zero is the binding-site start, but it need not coincide with
        # a geometric corner.  Rotating the closed perimeter places the whole duplex in
        # the middle of a straight edge and keeps it clear of both vertical sides.
        circle_points = perimeter_points[margin:] + perimeter_points[:margin]
        circle = FoldStrand(self.name, rotated, circle_points, self.circle_colour,
                            label="closed template")
        site_points = circle_points[:len(self.primer_sequence)]
        if any(abs(y - site_points[0][1]) > 1e-6 for _, y in site_points):
            raise ValueError("circular primer site must occupy one straight perimeter edge")
        # The rotated binding site is deliberately placed on the top edge.  A constant
        # inward offset makes an ordinary antiparallel duplex: parallel backbones and
        # vertical base-pair rungs, rather than a fan of misleading diagonal bonds.
        primer_points = tuple((x, y + 24.0) for x, y in reversed(site_points))
        primer = FoldStrand(self.primer_name, self.primer_sequence, tuple(primer_points),
                            self.primer_colour, label=self.primer_name)
        pairing = FoldPairing(circle.name, 0, primer.name, 0, len(primer.sequence),
                              f"{self.primer_name} binding site")
        closure = [(circle.name, len(circle.sequence) - 1, circle.name, 0,
                    "covalently closed circle")]
        return FoldLayout([circle, primer], [pairing], closure), width, height

    def folded_svg(self) -> str:
        layout, width, height = self.folded_layout()
        return layout.svg(width, height, self.name + " with bound primer")

    def folded_panel(self, caption: str = "") -> str:
        cap = f"<figcaption>{html.escape(caption)}</figcaption>" if caption else ""
        return (f'<figure class="chem-panel topology-panel">{cap}'
                f'<div class="diagram-scroll">{self.folded_svg()}</div></figure>')


@dataclass(frozen=True)
class PadlockCaptureTopology:
    """A molecular-inversion probe bound across one target gap.

    Probe order is 5′ ligation arm → unpaired backbone → extension arm 3′. Exact target
    sites are located from the arm sequences, so binding position and flexible backbone
    follow from sequence rather than caller-supplied coordinates.
    """

    name: str
    target_sequence: str
    ligation_arm: str
    backbone: str
    extension_arm: str

    def __post_init__(self) -> None:
        for attr in ("target_sequence", "ligation_arm", "backbone", "extension_arm"):
            object.__setattr__(self, attr, _clean(getattr(self, attr), attr))
        ext_site = self.target_sequence.find(revcomp(self.extension_arm))
        lig_site = self.target_sequence.find(revcomp(self.ligation_arm))
        if ext_site < 0 or self.target_sequence.find(revcomp(self.extension_arm), ext_site + 1) >= 0:
            raise ValueError("extension arm must have one exact target site")
        if lig_site < 0 or self.target_sequence.find(revcomp(self.ligation_arm), lig_site + 1) >= 0:
            raise ValueError("ligation arm must have one exact target site")
        ext_end = ext_site + len(self.extension_arm)
        lig_end = lig_site + len(self.ligation_arm)
        if not (ext_end <= lig_site or lig_end <= ext_site):
            raise ValueError("padlock extension and ligation arms must flank a non-overlapping gap")

    def folded_layout(self) -> FoldLayout:
        step = 12.0
        target_points = tuple((70 + step * i, 250.0)
                              for i in range(len(self.target_sequence)))
        target = FoldStrand("target", self.target_sequence, target_points, "var(--ink)",
                            label="target")
        ext_site = self.target_sequence.index(revcomp(self.extension_arm))
        lig_site = self.target_sequence.index(revcomp(self.ligation_arm))
        lig_points = tuple((target_points[lig_site + len(self.ligation_arm) - 1 - i][0], 226.0)
                           for i in range(len(self.ligation_arm)))
        ext_points = tuple((target_points[ext_site + len(self.extension_arm) - 1 - i][0], 226.0)
                           for i in range(len(self.extension_arm)))
        left, right = lig_points[-1], ext_points[0]
        extent = max(60.0, len(self.backbone) * 4.0)
        route = [left, (left[0], 176), ((left[0] + right[0]) / 2, 176 - extent / 4),
                 (right[0], 176), right]
        backbone_points = tuple(_sample_polyline(route, len(self.backbone)))
        probe_sequence = self.ligation_arm + self.backbone + self.extension_arm
        probe = FoldStrand("padlock probe", probe_sequence,
                           lig_points + backbone_points + ext_points, "var(--c-cbc)",
                           label="padlock probe")
        pairings = [
            FoldPairing(target.name, lig_site, probe.name, 0, len(self.ligation_arm),
                        "ligation-arm target site"),
            FoldPairing(target.name, ext_site, probe.name,
                        len(self.ligation_arm) + len(self.backbone),
                        len(self.extension_arm), "extension-arm target site"),
        ]
        straight = [StraightHint(target.name, 0, len(target.sequence) - 1)]
        return FoldLayout([target, probe], pairings, straight=straight)

    def folded_svg(self) -> str:
        width = max(620, 150 + len(self.target_sequence) * 12)
        return self.folded_layout().svg(width, 330, self.name + " gap capture")

    def folded_panel(self, caption: str = "") -> str:
        cap = f"<figcaption>{html.escape(caption)}</figcaption>" if caption else ""
        return (f'<figure class="chem-panel topology-panel">{cap}'
                f'<div class="diagram-scroll">{self.folded_svg()}</div></figure>')


@dataclass(frozen=True)
class AdjacentProbeDockingTopology:
    """Two adjacent target probes whose flexible tails form one docking duplex."""

    name: str
    left_target: str
    right_target: str
    left_spacer: str
    right_spacer: str
    left_tail: str
    right_tail: str
    docking_sequence: str

    def __post_init__(self) -> None:
        for attr in ("left_target", "right_target", "left_spacer", "right_spacer",
                     "left_tail", "right_tail", "docking_sequence"):
            object.__setattr__(self, attr, _clean(getattr(self, attr), attr))
        expected = revcomp(self.left_tail + self.right_tail)
        if self.docking_sequence != expected and not (
                set(self.docking_sequence) == {"N"} and len(self.docking_sequence) == len(expected)):
            raise ValueError("docking strand must complement both adjacent probe tails")

    def folded_layout(self) -> FoldLayout:
        step = 12.0
        target_sequence = self.left_target + self.right_target
        target_points = tuple((90 + step * i, 340.0) for i in range(len(target_sequence)))
        target = FoldStrand("target RNA", target_sequence, target_points, "var(--ink)",
                            "RNA", "target RNA")
        left_arm = tuple((target_points[len(self.left_target) - 1 - i][0], 316.0)
                         for i in range(len(self.left_target)))
        right_arm = tuple((target_points[len(target_sequence) - 1 - i][0], 316.0)
                          for i in range(len(self.right_target)))
        tail_total = len(self.left_tail) + len(self.right_tail)
        platform_x = 90 + (len(target_sequence) * step - tail_total * step) / 2
        left_tail_points = tuple((platform_x + step * i, 174.0)
                                 for i in range(len(self.left_tail)))
        right_tail_points = tuple((platform_x + step * (len(self.left_tail) + i), 174.0)
                                  for i in range(len(self.right_tail)))
        left_route = [left_arm[-1], (left_arm[-1][0] - 24, 275),
                      (left_tail_points[0][0] - 24, 214), left_tail_points[0]]
        right_route = [right_arm[-1], (right_arm[-1][0] - 24, 275),
                       (right_tail_points[0][0] - 24, 214), right_tail_points[0]]
        left_spacer_points = tuple(_sample_polyline(left_route, len(self.left_spacer)))
        right_spacer_points = tuple(_sample_polyline(right_route, len(self.right_spacer)))
        left_probe = FoldStrand("left target probe",
                                revcomp(self.left_target) + self.left_spacer + self.left_tail,
                                left_arm + left_spacer_points + left_tail_points,
                                "var(--c-r1)", label="left probe")
        right_probe = FoldStrand("right target probe",
                                 revcomp(self.right_target) + self.right_spacer + self.right_tail,
                                 right_arm + right_spacer_points + right_tail_points,
                                 "var(--c-r2)", label="right probe")
        combined_tail_points = left_tail_points + right_tail_points
        docking_points = tuple((x, 150.0) for x, _ in reversed(combined_tail_points))
        docking = FoldStrand("preamplifier docking strand", self.docking_sequence,
                             docking_points, "var(--c-tso)", label="preamplifier")
        left_tail_start = len(self.left_target) + len(self.left_spacer)
        right_tail_start = len(self.right_target) + len(self.right_spacer)
        pairings = [
            FoldPairing(target.name, 0, left_probe.name, 0, len(self.left_target),
                        "left target-probe site"),
            FoldPairing(target.name, len(self.left_target), right_probe.name, 0,
                        len(self.right_target), "right target-probe site"),
            FoldPairing(left_probe.name, left_tail_start, docking.name,
                        len(self.right_tail), len(self.left_tail),
                        "left half of docking site"),
            FoldPairing(right_probe.name, right_tail_start, docking.name, 0,
                        len(self.right_tail), "right half of docking site"),
        ]
        return FoldLayout([target, left_probe, right_probe, docking], pairings,
                          straight=[StraightHint(target.name, 0, len(target.sequence) - 1)])

    def folded_svg(self) -> str:
        width = max(760, 180 + (len(self.left_target) + len(self.right_target)) * 12)
        return self.folded_layout().svg(width, 420, self.name)

    def folded_panel(self, caption: str = "") -> str:
        cap = f"<figcaption>{html.escape(caption)}</figcaption>" if caption else ""
        return (f'<figure class="chem-panel topology-panel">{cap}'
                f'<div class="diagram-scroll">{self.folded_svg()}</div></figure>')


def _complements(first: str, second: str) -> bool:
    """Pairing equality in which N is an explicit source-unknown base."""
    if len(first) != len(second):
        return False
    return all(a == "N" or b == "N" or a == revcomp(b)
               for a, b in zip(first, reversed(second)))


@dataclass(frozen=True)
class PrimerTemplatedPadlockTopology:
    """Adjacent RNA probes with both padlock ends closed on the primer-probe tail."""

    name: str
    primer_target: str
    padlock_target: str
    target_gap: str
    primer_splint: str
    padlock_five_arm: str
    padlock_payload: str
    padlock_three_arm: str

    def __post_init__(self) -> None:
        for attr in ("primer_target", "padlock_target", "target_gap", "primer_splint",
                     "padlock_five_arm", "padlock_payload", "padlock_three_arm"):
            object.__setattr__(self, attr, _clean(getattr(self, attr), attr))
        if len(self.primer_splint) != (len(self.padlock_five_arm)
                                      + len(self.padlock_three_arm)):
            raise ValueError("primer splint must be exactly the two padlock docking arms")
        cut = len(self.padlock_five_arm)
        if not _complements(self.primer_splint[:cut], self.padlock_five_arm):
            raise ValueError("padlock 5-prime arm does not pair with the first splint half")
        if not _complements(self.primer_splint[cut:], self.padlock_three_arm):
            raise ValueError("padlock 3-prime arm does not pair with the second splint half")

    def folded_layout(self) -> FoldLayout:
        step = 12.0
        # Use the mirror-equivalent planar orientation with the padlock-binding
        # RNA site at the left.  Mirroring changes only page coordinates, while
        # preserving strand order and every molecular adjacency.
        target_sequence = self.padlock_target + self.target_gap + self.primer_target
        target_points = tuple((90 + step * i, 390.0) for i in range(len(target_sequence)))
        target = FoldStrand("target RNA", target_sequence, target_points, "var(--ink)",
                            "RNA", "target RNA")

        primer_site = len(self.padlock_target) + len(self.target_gap)
        primer_arm = tuple((target_points[primer_site + len(self.primer_target) - 1 - i][0],
                            366.0)
                           for i in range(len(self.primer_target)))
        platform_x = primer_arm[-1][0]
        splint_points = tuple((platform_x + step * i, 170.0)
                              for i in range(len(self.primer_splint)))
        primer = FoldStrand(
            "primer probe", revcomp(self.primer_target) + self.primer_splint,
            primer_arm + splint_points, "var(--c-r1)", label="primer probe / RCA primer")

        cut = len(self.padlock_five_arm)
        five_points = tuple((splint_points[cut - 1 - i][0], 146.0)
                            for i in range(len(self.padlock_five_arm)))
        padlock_site = 0
        target_arm = tuple((target_points[len(self.padlock_target) - 1 - i][0],
                            366.0) for i in range(len(self.padlock_target)))
        three_points = tuple((splint_points[-1 - i][0], 146.0)
                             for i in range(len(self.padlock_three_arm)))
        route = [target_arm[-1], (40.0, 300.0), (40.0, 108.0),
                 (splint_points[-1][0] + 68.0, 108.0),
                 (splint_points[-1][0] + 68.0, 132.0), three_points[0]]
        payload_points = tuple(_sample_polyline(route, len(self.padlock_payload)))
        padlock = FoldStrand(
            "padlock probe",
            self.padlock_five_arm + revcomp(self.padlock_target)
            + self.padlock_payload + self.padlock_three_arm,
            five_points + target_arm + payload_points + three_points,
            "var(--c-cbc)", label="5′-phosphorylated padlock")

        primer_splint_start = len(self.primer_target)
        padlock_target_start = len(self.padlock_five_arm)
        padlock_three_start = (padlock_target_start + len(self.padlock_target)
                               + len(self.padlock_payload))
        pairings = [
            FoldPairing(target.name, primer_site, primer.name, 0, len(self.primer_target),
                        "primer-probe target site"),
            FoldPairing(target.name, padlock_site, padlock.name, padlock_target_start,
                        len(self.padlock_target), "padlock target site"),
            FoldPairing(primer.name, primer_splint_start, padlock.name, 0,
                        len(self.padlock_five_arm), "padlock 5-prime docking arm"),
            FoldPairing(primer.name, primer_splint_start + cut, padlock.name,
                        padlock_three_start, len(self.padlock_three_arm),
                        "padlock 3-prime docking arm"),
        ]
        closure = [(padlock.name, len(padlock.sequence) - 1, padlock.name, 0,
                    "T4 DNA ligase padlock closure")]
        return FoldLayout([target, primer, padlock], pairings, closure,
                          [StraightHint(target.name, 0, len(target.sequence) - 1)])

    def folded_svg(self) -> str:
        width = max(900, 220 + (len(self.primer_target) + len(self.target_gap)
                               + len(self.padlock_target)) * 12)
        return self.folded_layout().svg(width, 460, self.name)

    def folded_panel(self, caption: str = "") -> str:
        cap = f"<figcaption>{html.escape(caption)}</figcaption>" if caption else ""
        return (f'<figure class="chem-panel topology-panel">{cap}'
                f'<div class="diagram-scroll">{self.folded_svg()}</div></figure>')


@dataclass(frozen=True)
class SplintProductSegment:
    """One named region in a product assembled on an interrupted splint."""

    name: str
    sequence: str
    docking: bool = False
    junction_before: bool = False


@dataclass(frozen=True)
class InterruptedSplintTopology:
    """A ligated product aligned across nucleotide chunks separated by non-base spacers."""

    name: str
    splint_chunks: tuple[tuple[str, str], ...]
    product_segments: tuple[SplintProductSegment, ...]
    spacer_label: str = "C3"

    def __post_init__(self) -> None:
        if len(self.splint_chunks) < 2:
            raise ValueError("an interrupted splint needs at least two nucleotide chunks")
        clean_chunks = tuple((name, _clean(sequence, name))
                             for name, sequence in self.splint_chunks)
        clean_segments = tuple(SplintProductSegment(
            segment.name, _clean(segment.sequence, segment.name),
            segment.docking, segment.junction_before)
            for segment in self.product_segments)
        object.__setattr__(self, "splint_chunks", clean_chunks)
        object.__setattr__(self, "product_segments", clean_segments)
        splint = "".join(sequence for _, sequence in clean_chunks)
        previous = len(splint)
        for segment in clean_segments:
            if not segment.docking:
                continue
            site = splint.find(revcomp(segment.sequence))
            if site < 0 or splint.find(revcomp(segment.sequence), site + 1) >= 0:
                raise ValueError(f"{segment.name} must have one exact splint site")
            if site + len(segment.sequence) > previous:
                raise ValueError("product docking regions must traverse the splint antiparallel")
            previous = site

    def folded_layout(self) -> tuple[FoldLayout, tuple[Point, ...], tuple[Point, ...], int]:
        step, gap, x0, y = 12.0, 42.0, 70.0, 300.0
        splint_sequence = "".join(sequence for _, sequence in self.splint_chunks)
        splint_points: list[Point] = []
        spacer_points = []
        x = x0
        for chunk_index, (_, sequence) in enumerate(self.splint_chunks):
            for _ in sequence:
                splint_points.append((x, y)); x += step
            if chunk_index + 1 < len(self.splint_chunks):
                spacer_points.append((x - step / 2 + gap / 2, y))
                x += gap
        splint = FoldStrand("interrupted splint", splint_sequence, tuple(splint_points),
                            "var(--c-s5)", label="Splint6")

        product_sequence = "".join(segment.sequence for segment in self.product_segments)
        product_points: list[Point | None] = [None] * len(product_sequence)
        junction_indices = []
        cursor = 0
        for segment in self.product_segments:
            if segment.junction_before:
                if cursor == 0:
                    raise ValueError("a ligation junction cannot precede the product")
                junction_indices.append(cursor)
            if segment.docking:
                site = splint_sequence.index(revcomp(segment.sequence))
                for offset in range(len(segment.sequence)):
                    sx, _ = splint_points[site + len(segment.sequence) - 1 - offset]
                    product_points[cursor + offset] = (sx, y - 26)
            cursor += len(segment.sequence)

        index = 0
        while index < len(product_points):
            if product_points[index] is not None:
                index += 1; continue
            start = index
            while index < len(product_points) and product_points[index] is None:
                index += 1
            count = index - start
            before = product_points[start - 1] if start else None
            after = product_points[index] if index < len(product_points) else None
            if before and after:
                middle = ((before[0] + after[0]) / 2, y - 92)
                routed = _sample_polyline([before, middle, after], count)
            elif after:
                routed = tuple((after[0] + step * (count - offset), y - 70)
                               for offset in range(count))
            elif before:
                routed = tuple((before[0] - step * (offset + 1), y - 70)
                               for offset in range(count))
            else:
                raise ValueError("the product must contain at least one docking region")
            product_points[start:index] = routed

        product = FoldStrand("assembled product", product_sequence,
                             tuple(product_points), "var(--c-cbc)", label="ligated product")
        pairings = []
        product_cursor = 0
        for segment in self.product_segments:
            if segment.docking:
                site = splint_sequence.index(revcomp(segment.sequence))
                pairings.append(FoldPairing(splint.name, site, product.name, product_cursor,
                                            len(segment.sequence), segment.name + " docking"))
            product_cursor += len(segment.sequence)
        layout = FoldLayout([splint, product], pairings)
        junction_points = tuple(
            ((product.points[i - 1][0] + product.points[i][0]) / 2,
             (product.points[i - 1][1] + product.points[i][1]) / 2)
            for i in junction_indices)
        return layout, tuple(spacer_points), junction_points, int(x + 60)

    def folded_svg(self) -> str:
        layout, spacers, junctions, width = self.folded_layout()
        svg = layout.svg(max(900, width), 370, self.name)
        additions = []
        for x, y in spacers:
            additions.append(f'<text x="{x:.1f}" y="{y + 5:.1f}" text-anchor="middle" '
                             f'class="topology-label">{html.escape(self.spacer_label)}</text>')
        for x, y in junctions:
            additions.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="5" '
                             'class="topology-junction"><title>T4 DNA ligase junction</title></circle>')
        return svg.replace("</svg>", "".join(additions) + "</svg>")

    def folded_panel(self, caption: str = "") -> str:
        cap = f"<figcaption>{html.escape(caption)}</figcaption>" if caption else ""
        return (f'<figure class="chem-panel topology-panel">{cap}'
                f'<div class="diagram-scroll">{self.folded_svg()}</div></figure>')
