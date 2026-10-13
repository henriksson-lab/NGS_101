"""Reusable molecular states for isothermal amplification chemistries.

These helpers describe topology, not a fictitious single consensus sequence.  LAMP and
multiple-displacement products are populations of branched molecules, so callers supply
the sequence-bearing primers separately and use these states for the reaction workflow.
"""
from __future__ import annotations

from dataclasses import dataclass
import html

from chemdraw import Scene
from protocol_model import molecule, role


@dataclass(frozen=True)
class LampTopology:
    """Deterministic role-level drawings of the two characteristic LAMP topologies.

    A LAMP reaction is a heterogeneous population, not one reconstructable sequence.
    These views therefore draw connectivity and repeated molecular roles without
    assigning bases or pretending that branch count and placement are universal.
    """

    target: str = "target-derived DNA"
    barcode: str = "FIP sample barcode"

    @staticmethod
    def _panel(svg: str, caption: str) -> str:
        return (f'<figure class="chem-panel topology-panel"><figcaption>'
                f'{html.escape(caption)}</figcaption><div class="diagram-scroll">'
                f'{svg}</div></figure>')

    @staticmethod
    def _label(x: float, y: float, text: str, *, cls: str = "topology-label") -> str:
        return (f'<text x="{x}" y="{y}" text-anchor="middle" class="{cls}">'
                f'{html.escape(text)}</text>')

    def stem_loop_svg(self) -> str:
        """Draw one strand folded into the first double-ended stem-loop precursor."""
        target = html.escape(self.target)
        parts = [
            '<svg xmlns="http://www.w3.org/2000/svg" class="chem-svg topology-svg" '
            'width="940" height="270" viewBox="0 0 940 270" role="img" '
            'aria-label="LAMP double-ended stem-loop intermediate; role-level topology">',
            '<title>LAMP stem-loop intermediate; nucleotide order varies with target and primer set</title>',
            # One continuous strand, opened only at the adjacent 5-prime/3-prime ends.
            '<path d="M480,88 L720,88 C790,88 790,182 720,182 '
            'L180,182 C110,182 110,88 180,88 L450,88" fill="none" '
            'stroke="var(--c-r1)" stroke-width="5" stroke-linecap="round"/>',
            '<g stroke="var(--note-rule)" stroke-width="1.5">',
            *[f'<line x1="{x}" y1="94" x2="{x}" y2="176"/>'
              for x in range(190, 721, 38)],
            '</g>',
            '<circle cx="480" cy="88" r="4" fill="var(--c-r1)"/>',
            '<circle cx="450" cy="88" r="4" fill="var(--c-r1)"/>',
            self._label(486, 70, "5′"), self._label(444, 70, "3′"),
            self._label(145, 140, "FIP-derived"), self._label(145, 158, "terminal loop"),
            self._label(755, 140, "BIP-derived"), self._label(755, 158, "terminal loop"),
            self._label(450, 128, target),
            self._label(450, 148, "paired inverted target copies", cls="topology-comment"),
            self._label(450, 230, "symbolic topology — bases and loop lengths depend on the assay",
                        cls="topology-comment"),
            '</svg>',
        ]
        return "".join(parts)

    def branched_svg(self) -> str:
        """Draw a representative non-crossing branched concatemer population member."""
        target = html.escape(self.target)
        barcode = html.escape(self.barcode)
        parts = [
            '<svg xmlns="http://www.w3.org/2000/svg" class="chem-svg topology-svg" '
            'width="1040" height="360" viewBox="0 0 1040 360" role="img" '
            'aria-label="branched LAMP concatemer; role-level topology">',
            '<title>Representative branched LAMP product; branch number and order vary</title>',
            '<path d="M55,180 L985,180" fill="none" stroke="var(--c-r1)" '
            'stroke-width="5" stroke-linecap="round"/>',
            # Each branch leaves and rejoins the main strand locally. Alternating lanes
            # make the fixed embedding planar without implying a biological orientation.
            '<path d="M236,180 L236,98 C236,48 304,48 304,98 L304,180" '
            'fill="none" stroke="var(--c-r2)" stroke-width="5"/>',
            '<path d="M486,180 L486,262 C486,312 554,312 554,262 L554,180" '
            'fill="none" stroke="var(--c-r2)" stroke-width="5"/>',
            '<path d="M736,180 L736,98 C736,48 804,48 804,98 L804,180" '
            'fill="none" stroke="var(--c-r2)" stroke-width="5"/>',
            '<g stroke="var(--note-rule)" stroke-width="1.5">',
            *[f'<line x1="{a}" y1="{y}" x2="{b}" y2="{y}"/>'
              for a, b, ys in ((242, 298, (116, 140, 164)),
                                (492, 548, (200, 224, 248)),
                                (742, 798, (116, 140, 164)))
              for y in ys],
            '</g>',
        ]
        for centre in (145, 395, 645, 895):
            parts.extend([
                f'<rect x="{centre-82}" y="158" width="112" height="44" rx="8" '
                'fill="var(--surface)" stroke="var(--c-me)" stroke-width="2">'
                f'<title>{target}; exact bases depend on target</title></rect>',
                self._label(centre - 26, 176, "target-derived"),
                self._label(centre - 26, 193, "repeat"),
                f'<rect x="{centre+34}" y="158" width="48" height="44" rx="8" '
                'fill="var(--surface)" stroke="var(--c-cbc)" stroke-width="2">'
                f'<title>{barcode}; bases may be proprietary or assay-specific</title></rect>',
                self._label(centre + 58, 184, "BC"),
            ])
        parts.extend([
            self._label(270, 34, "stem-loop branch"),
            self._label(520, 338, "displaced stem-loop branch"),
            self._label(770, 34, "stem-loop branch"),
            self._label(520, 228, "repeated inverted units along a displaced concatemer",
                        cls="topology-comment"),
            self._label(520, 252, "representative population member; branch count and positions vary",
                        cls="topology-comment"),
            '</svg>',
        ])
        return "".join(parts)

    def panels(self) -> tuple[tuple[str, str], ...]:
        return (
            ("Stem-loop precursor", self._panel(
                self.stem_loop_svg(),
                "The first displaced LAMP product folds its inverted primer-derived ends "
                "into a double-ended stem-loop. This is molecular connectivity, not a "
                "single claimed nucleotide sequence.")),
            ("Branched LAMP product", self._panel(
                self.branched_svg(),
                "Continued strand displacement produces a population of concatemers with "
                "repeated target/barcode units and stem-loop branches. One planar member is "
                "shown; branch count, order and exact bases are deliberately not asserted.")),
        )


def lamp_topology(target: str = "target-derived DNA",
                  barcode: str = "FIP sample barcode") -> LampTopology:
    """Return the shared symbolic topology used by LAMP sequencing protocols."""
    if not target.strip() or not barcode.strip():
        raise ValueError("LAMP topology roles must be named")
    return LampTopology(target, barcode)


def lamp_stem_loop(target: str = "target-derived DNA") -> Scene:
    """First self-priming dumbbell/stem-loop intermediate produced by inner primers."""
    return molecule("LAMP stem-loop intermediate", (
        role("FIP-derived inverted repeat"), role(target),
        role("BIP-derived inverted repeat"), role("self-paired terminal loops")))


def lamp_product(target: str = "target-derived DNA",
                 barcode: str = "FIP sample barcode") -> Scene:
    """Branched concatemer population after loop-primer-assisted displacement."""
    return molecule("branched LAMP concatemer population", (
        role("repeated inverted target units"), role(barcode, tag="cbc"),
        role("stem loops and displaced branches"), role(target)))


def mda_product(label: str = "branched multiple-displacement product") -> Scene:
    """Random-primer/phi29 product without implying one linear final molecule."""
    return molecule(label, (
        role("long copied genomic tracts"), role("displaced single strands"),
        role("new random-primer extension branches")))


def spia_product() -> Scene:
    """Ribo-SPIA cDNA population made through repeated RNase-H priming cycles."""
    return molecule("Ribo-SPIA amplified cDNA population", (
        role("repeatedly displaced transcript-derived cDNA copies"),
        role("complementary DNA products"), role("SPIA handle remnant")))


def proximity_extension_pair() -> Scene:
    """Two antibody probes co-bound to one analyte before polymerase extension."""
    return molecule("antigen-bound PEA probe pair", (
        role("forward antibody + oligo"), role("target protein"),
        role("reverse antibody + oligo"), role("paired hybridization arms")))


def proximity_extension_product() -> Scene:
    """Reporter created when one probe arm templates extension of its neighbour."""
    return molecule("proximity-extension reporter", (
        role("P5 and Read 1 precursor"), role("forward assay barcode", tag="cbc"),
        role("copied hybridization junction"), role("reverse assay barcode", tag="cbc"),
        role("common reverse-PCR handle")))
