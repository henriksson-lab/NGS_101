"""Deterministic role-level topology for proximity assays with undisclosed oligos.

These drawings encode molecular connectivity, not thermodynamic folds or guessed bases.
Protocol modules may supply exact published binding regions; everything else remains a
named role.  Fixed planar embeddings make antibody, analyte, probe and reporter ownership
clear without crossings.
"""
from __future__ import annotations

from dataclasses import dataclass
import html
import math

from chemdraw import Segment


def _label(x: float, y: float, text: str, *, cls: str = "topology-label",
           anchor: str = "middle") -> str:
    return (f'<text x="{x}" y="{y}" text-anchor="{anchor}" class="{cls}">'
            f'{html.escape(text)}</text>')


def _panel(svg: str, caption: str) -> str:
    return (f'<figure class="chem-panel topology-panel"><figcaption>'
            f'{html.escape(caption)}</figcaption><div class="diagram-scroll">'
            f'{svg}</div></figure>')


def _antibody(x: float, y: float, label: str, colour: str) -> str:
    """One Y-shaped antibody with its oligo conjugated at the Fc stem."""
    return "".join((
        f'<g><title>{html.escape(label)} with covalently attached DNA probe</title>',
        f'<path d="M{x},{y+68} L{x},{y+28} M{x},{y+28} L{x-28},{y} '
        f'M{x},{y+28} L{x+28},{y}" fill="none" stroke="{colour}" '
        'stroke-width="7" stroke-linecap="round"/>',
        f'<circle cx="{x}" cy="{y+68}" r="5" fill="{colour}"/>',
        '</g>',
        _label(x, y - 12, label),
    ))


@dataclass(frozen=True)
class ProximityExtensionTopology:
    """Two antibody probes co-bind one analyte and create an extension reporter."""

    analyte: str = "target protein"
    left_probe: str = "forward antibody–oligo"
    right_probe: str = "reverse antibody–oligo"

    def __post_init__(self) -> None:
        if not all(x.strip() for x in (self.analyte, self.left_probe, self.right_probe)):
            raise ValueError("proximity-extension topology roles must be named")

    def svg(self) -> str:
        parts = [
            '<svg xmlns="http://www.w3.org/2000/svg" class="chem-svg topology-svg" '
            'width="1080" height="370" viewBox="0 0 1080 370" role="img" '
            'aria-label="antibody proximity-extension topology">',
            '<title>Matched antibody oligos co-bind one analyte, hybridize and template extension</title>',
            _label(260, 28, "1. Matched probes co-bind one analyte", cls="topology-heading"),
            _label(810, 28, "2. Proximal arms hybridize and extend", cls="topology-heading"),
        ]
        for offset, paired in ((0, False), (550, True)):
            cx = offset + 260
            parts.extend([
                f'<ellipse cx="{cx}" cy="82" rx="72" ry="30" fill="var(--surface)" '
                'stroke="var(--note-rule)" stroke-width="3">'
                f'<title>{html.escape(self.analyte)}</title></ellipse>',
                _label(cx, 87, self.analyte),
                _antibody(cx - 92, 112, self.left_probe, "var(--c-r1)"),
                _antibody(cx + 92, 112, self.right_probe, "var(--c-r2)"),
                f'<path d="M{cx-92},208 L{cx-92},252 L{cx-20},252" fill="none" '
                'stroke="var(--c-r1)" stroke-width="5" stroke-linecap="round"/>',
                f'<path d="M{cx+92},208 L{cx+92},280 L{cx+20},280" fill="none" '
                'stroke="var(--c-r2)" stroke-width="5" stroke-linecap="round"/>',
                _label(cx - 105, 238, "forward barcode", anchor="end"),
                _label(cx + 105, 303, "reverse barcode", anchor="start"),
            ])
            if paired:
                parts.append('<g stroke="var(--note-rule)" stroke-width="1.5">')
                parts.extend(f'<line x1="{x}" y1="257" x2="{x}" y2="275"/>'
                             for x in range(int(cx - 14), int(cx + 15), 7))
                parts.extend([
                    '</g>',
                    f'<path d="M{cx-20},252 L{cx+72},252" fill="none" '
                    'stroke="var(--c-cbc)" stroke-width="8" stroke-linecap="round"/>',
                    f'<path d="M{cx+62},246 l12,6 -12,6 z" fill="var(--c-cbc)"/>',
                    _label(cx + 12, 232, "polymerase extension"),
                    _label(cx, 332, "paired assay barcodes become one amplifiable reporter",
                           cls="topology-comment"),
                ])
            else:
                parts.append(_label(cx, 332, "unbound or mismatched probes do not form the paired template",
                                    cls="topology-comment"))
        parts.append('</svg>')
        return "".join(parts)

    def panels(self) -> tuple[tuple[str, str], ...]:
        return (("Proximity-extension topology", _panel(
            self.svg(),
            "Two antibody–oligo conjugates must bind the same analyte closely enough for "
            "their terminal arms to pair. Polymerase then copies across the paired arms, "
            "joining the two assay identities in one reporter. Probe bases are proprietary "
            "and are therefore shown by role only.")),)


@dataclass(frozen=True)
class ProximityLigationTopology:
    """Two antibody probes are joined only through a shared splint/connector."""

    left_barcode: str = "antibody A barcode"
    right_barcode: str = "antibody B barcode"
    connector: str = "universal connector"

    def __post_init__(self) -> None:
        if not all(x.strip() for x in (self.left_barcode, self.right_barcode,
                                       self.connector)):
            raise ValueError("proximity-ligation topology roles must be named")

    def svg(self) -> str:
        parts = [
            '<svg xmlns="http://www.w3.org/2000/svg" class="chem-svg topology-svg" '
            'width="940" height="390" viewBox="0 0 940 390" role="img" '
            'aria-label="antibody proximity ligation through a universal connector">',
            '<title>Adjacent antibody probes hybridize to a common connector and are ligated</title>',
            '<ellipse cx="375" cy="72" rx="55" ry="27" fill="var(--surface)" '
            'stroke="var(--note-rule)" stroke-width="3"/>',
            '<ellipse cx="565" cy="72" rx="55" ry="27" fill="var(--surface)" '
            'stroke="var(--note-rule)" stroke-width="3"/>',
            _label(375, 77, "protein A"), _label(565, 77, "protein B"),
            _antibody(390, 112, "probe A antibody", "var(--c-r1)"),
            _antibody(550, 112, "probe B antibody", "var(--c-r2)"),
            '<path d="M390,208 L390,268 L470,268" fill="none" '
            'stroke="var(--c-r1)" stroke-width="6" stroke-linecap="round"/>',
            '<path d="M550,208 L550,268 L470,268" fill="none" '
            'stroke="var(--c-r2)" stroke-width="6" stroke-linecap="round"/>',
            '<path d="M430,296 L510,296" fill="none" stroke="var(--c-tso)" '
            'stroke-width="8" stroke-linecap="round"/>',
            '<g stroke="var(--note-rule)" stroke-width="1.5">',
            *[f'<line x1="{x}" y1="274" x2="{x}" y2="290"/>'
              for x in range(434, 511, 8)],
            '</g>',
            '<circle cx="470" cy="268" r="6" class="topology-junction">'
            '<title>connector-templated ligation junction</title></circle>',
            _label(350, 250, self.left_barcode, anchor="end"),
            _label(590, 250, self.right_barcode, anchor="start"),
            _label(470, 323, self.connector),
            _label(470, 355, "the connector splints one A–B nick; it is not inserted into the reporter",
                   cls="topology-comment"),
            '</svg>',
        ]
        return "".join(parts)

    def panels(self) -> tuple[tuple[str, str], ...]:
        return (("Connector-templated proximity ligation", _panel(
            self.svg(),
            "A common connector hybridizes across the adjacent ends of two proximal "
            "antibody probes and splints their single ligation junction. The covalent "
            "reporter comprises probe A joined directly to probe B; the connector is "
            "shown as a hybridized structural reagent, not inserted sequence. "
            "Undisclosed bases are not invented.")),)


@dataclass(frozen=True)
class ProximityNetworkTopology:
    """Two rolling-circle nodes bridged into one type-1/type-2 proximity product."""

    type1_binding_site: str
    type2_binding_site: str
    event_identifier_length: int

    def __post_init__(self) -> None:
        for name, sequence in (("type-1 linker-binding site", self.type1_binding_site),
                               ("type-2 linker-binding site", self.type2_binding_site)):
            Segment(name, sequence)  # shared sequence validation and hover semantics
        if self.event_identifier_length < 1:
            raise ValueError("a proximity link needs a positive event-identifier length")

    @staticmethod
    def _spiral(cx: float, cy: float, phase: float) -> str:
        points = []
        for i in range(80):
            angle = phase + i * math.pi / 10
            radius = 8 + i * 0.9
            points.append((cx + radius * math.cos(angle), cy + radius * math.sin(angle)))
        return " ".join(f"{x:.1f},{y:.1f}" for x, y in points)

    def svg(self) -> str:
        site1 = html.escape(Segment("type-1 linker-binding site",
                                    self.type1_binding_site).hover_text())
        site2 = html.escape(Segment("type-2 linker-binding site",
                                    self.type2_binding_site).hover_text())
        n1, n2 = len(self.type1_binding_site), len(self.type2_binding_site)
        parts = [
            '<svg xmlns="http://www.w3.org/2000/svg" class="chem-svg topology-svg" '
            'width="1120" height="430" viewBox="0 0 1120 430" role="img" '
            'aria-label="Pixelgen proximity link between two rolling-circle nodes">',
            '<title>Type-1 and type-2 rolling-circle nodes bridged by one proximity linker</title>',
            f'<polyline points="{self._spiral(235, 182, 0)}" fill="none" '
            'stroke="var(--c-r1)" stroke-width="5" stroke-linecap="round"/>',
            f'<polyline points="{self._spiral(885, 182, math.pi)}" fill="none" '
            'stroke="var(--c-r2)" stroke-width="5" stroke-linecap="round"/>',
            _label(235, 42, "type-1 RCA node", cls="topology-heading"),
            _label(885, 42, "type-2 RCA node", cls="topology-heading"),
            _label(235, 190, "PID-1 + UMI-1 repeats"),
            _label(885, 190, "PID-2 + UMI-2 repeats"),
            '<path d="M315,182 L370,182" stroke="var(--c-r1)" stroke-width="5"/>',
            '<path d="M805,182 L750,182" stroke="var(--c-r2)" stroke-width="5"/>',
            '<g><title>' + site1 + '</title><rect x="370" y="158" width="190" height="48" '
            'rx="8" fill="var(--surface)" stroke="var(--c-me)" stroke-width="2"/></g>',
            f'<text x="465" y="187" text-anchor="middle" class="fold-base">'
            f'{html.escape(self.type1_binding_site)}</text>',
            f'<g><title>unique event identifier — {self.event_identifier_length} nt</title>'
            '<rect x="560" y="158" width="90" height="48" rx="8" fill="var(--surface)" '
            'stroke="var(--c-cbc)" stroke-width="2"/></g>',
            _label(605, 181, f"UEI {self.event_identifier_length}N"),
            '<g><title>' + site2 + '</title><rect x="650" y="158" width="100" height="48" '
            'rx="8" fill="var(--surface)" stroke="var(--c-me)" stroke-width="2"/></g>',
            f'<text x="700" y="187" text-anchor="middle" class="fold-base">'
            f'{html.escape(self.type2_binding_site)}</text>',
            '<circle cx="370" cy="182" r="6" class="topology-junction">'
            '<title>type-1 gap-fill/ligation boundary</title></circle>',
            '<circle cx="750" cy="182" r="6" class="topology-junction">'
            '<title>type-2 gap-fill/ligation boundary</title></circle>',
            _label(465, 226, f"type-1 binding site ({n1} nt)"),
            _label(700, 226, f"type-2 binding site ({n2} nt)"),
            _label(560, 282, "gap fill and ligation retain one type-1/type-2 neighborhood edge",
                   cls="topology-comment"),
            _label(560, 310, "RCA coils are symbolic population nodes; their repeat count is not fixed",
                   cls="topology-comment"),
            '<path d="M235,300 L235,352" stroke="var(--ink-muted)" stroke-width="3"/>',
            '<path d="M885,300 L885,352" stroke="var(--ink-muted)" stroke-width="3"/>',
            _label(235, 376, "antibody-bound type-1 padlock origin"),
            _label(885, 376, "antibody-bound type-2 padlock origin"),
            '</svg>',
        ]
        return "".join(parts)

    def panels(self) -> tuple[tuple[str, str], ...]:
        return (("RCA-node proximity-link topology", _panel(
            self.svg(),
            "Each antibody-bound padlock becomes a localized rolling-circle node. A "
            "linker can bridge a neighboring type-1 and type-2 node; gap fill and ligation "
            "produce one edge carrying both PID–UMI identities and a unique event identifier. "
            "The two displayed linker-binding sequences are source-exact; undisclosed "
            "padlock and antibody-conjugate bases remain symbolic.")),)
