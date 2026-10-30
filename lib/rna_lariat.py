"""Declared, nucleotide-resolved RNA lariat connectivity.

The caller supplies RNA and a zero-based branch coordinate. This describes an
unmodified RNA backbone with one additional 2′–5′ bond; it does not predict a
splice site, splicing, nuclease resistance, reverse transcription or folding.
"""
from __future__ import annotations

from dataclasses import dataclass
import html
import math

from topology import FoldLayout, FoldStrand
from chemdraw import routed_strand_polygon


@dataclass(frozen=True)
class RNAOxygenBond:
    phosphate_nt: int
    oxygen_nt: int
    oxygen_position: int


@dataclass(frozen=True)
class RNALariatTopology:
    sequence: str
    branch_index: int
    name: str = 'RNA lariat'

    def __post_init__(self):
        if not isinstance(self.sequence, str):
            raise ValueError('RNA sequence must be text')
        object.__setattr__(self, 'sequence', ''.join(self.sequence.split()).upper())
        self.validate_connectivity()

    @property
    def bonds(self) -> tuple[RNAOxygenBond, ...]:
        return (tuple(RNAOxygenBond(i + 1, i, 3) for i in range(len(self.sequence) - 1))
                + (RNAOxygenBond(0, self.branch_index, 2),))

    def validate_connectivity(self) -> None:
        if len(self.sequence) < 2 or set(self.sequence) - set('ACGUN'):
            raise ValueError('lariat reference needs unmodified RNA letters A/C/G/U/N')
        if type(self.branch_index) is not int or not 0 < self.branch_index < len(self.sequence):
            raise ValueError('branch must be a distinct nucleotide inside the supplied RNA')
        expected = tuple(RNAOxygenBond(i + 1, i, 3) for i in range(len(self.sequence) - 1))
        expected += (RNAOxygenBond(0, self.branch_index, 2),)
        if self.bonds != expected:
            raise ValueError('backbone must retain every ordinary bond and the declared 2′–5′ link')
        # Different ribose ports are real endpoints: the branch's 3′ backbone
        # and additional 2′ bond must never be conflated into one attachment.
        if len({b.phosphate_nt for b in self.bonds}) != len(self.bonds):
            raise ValueError('a phosphate has more than one declared backbone attachment')
        if len({(b.oxygen_nt, b.oxygen_position) for b in self.bonds}) != len(self.bonds):
            raise ValueError('an oxygen has more than one declared backbone attachment')

    @property
    def loop_sequence(self) -> str:
        return self.sequence[:self.branch_index + 1]

    @property
    def tail_sequence(self) -> str:
        """Ordinary backbone downstream of the branch; excludes the branch base."""
        return self.sequence[self.branch_index + 1:]

    @property
    def five_prime_is_free(self) -> bool:
        return not any(b.phosphate_nt == 0 for b in self.bonds)

    @property
    def free_three_prime_index(self) -> int | None:
        last = len(self.sequence) - 1
        return None if any(b.oxygen_nt == last and b.oxygen_position == 3 for b in self.bonds) else last

    def layout(self) -> FoldLayout:
        self.validate_connectivity()
        loop_n = self.branch_index + 1
        radius = max(70.0, loop_n * 16 / (2 * math.pi))
        cx, cy = radius + 22, radius + 60
        points = []
        for i in range(loop_n):
            angle = 2 * math.pi * (i + 1) / loop_n
            points.append((cx + radius * math.cos(angle), cy + radius * math.sin(angle)))
        for i in range(len(self.tail_sequence)):
            points.append((cx + radius + (i + 1) * 22, cy))
        # Reuse the shared topology geometry and its crossing validator. The
        # generic renderer labels free termini, so lariat termini render below.
        return FoldLayout([FoldStrand(self.name, self.sequence, tuple(points), polymer='RNA')], [],
                          extra_bonds=[(self.name, self.branch_index, self.name, 0, '2′–5′ phosphodiester')])

    def svg(self, *, caption: str) -> str:
        layout = self.layout()
        points = layout.strands[self.name].points
        width = max(380, math.ceil(max(x for x, _ in points) + 50))
        height = math.ceil(max(y for _, y in points) + 55)
        e = html.escape
        out = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
               f'width="{width}" height="{height}" role="img" aria-label="{e(caption, quote=True)}">',
               f'<title>{e(caption)}</title>',
               '<style>text{font-family:system-ui,sans-serif;font-size:12px;fill:#24394b} '
               '.nt{font-family:monospace;font-size:12px}.bond-label{fill:#a33236}</style>',
               f'<text x="16" y="24">{e(caption)}</text>']
        shape = routed_strand_polygon(points, width=7.0, arrow_at='end')
        out.append(f'<polygon points="{shape}" fill="#edf3f8" stroke="#47677f" '
                   'stroke-width="1.4" stroke-linejoin="round" data-backbone="3-prime-5-prime"/>')
        bx, by = points[self.branch_index]
        sx, sy = points[0]
        out.append('<g><title>Branch 2′ oxygen to splice-site 5′ phosphate</title>'
                   f'<path d="M{bx:.2f},{by:.2f} L{sx:.2f},{sy:.2f}" '
                   'fill="none" stroke="#a33236" stroke-width="3" data-oxygen-position="2" '
                   'data-phosphate-nt="0"/></g>')
        for i, (base, (x, y)) in enumerate(zip(self.sequence, points)):
            out.append(f'<g data-nt="{i}"><title>{base}, RNA nucleotide {i+1}</title>'
                       f'<text class="nt" x="{x:.2f}" y="{y+4:.2f}" text-anchor="middle">{base}</text></g>')
        out.append(f'<text class="bond-label" x="{bx+8:.2f}" y="{by+35:.2f}">2′–5′ bond; branch nt {self.branch_index+1}</text>')
        out.append(f'<text x="{sx-48:.2f}" y="{sy+26:.2f}">5′ splice site, joined</text>')
        x3, y3 = points[self.free_three_prime_index]
        out.append(f'<text x="{x3+8:.2f}" y="{y3-13:.2f}" data-end="3-prime-free">3′ OH</text>')
        out.append('</svg>')
        return ''.join(out)
