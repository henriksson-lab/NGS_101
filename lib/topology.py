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

from chemdraw import Segment, revcomp


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

    def panel(self, caption: str = "") -> str:
        cap = f"<figcaption>{html.escape(caption)}</figcaption>" if caption else ""
        return (f'<figure class="chem-panel topology-panel">{cap}'
                f'<div class="diagram-scroll">{self.svg()}</div></figure>')
