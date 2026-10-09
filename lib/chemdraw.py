"""
chemdraw -- build character-aligned DNA construct diagrams in the scg_lib_structs idiom.

The point of this module is correctness. Hand-aligning a dozen near-identical duplex
ladders is exactly the kind of work that goes quietly wrong: a reverse complement written
backwards, a primer arrow two columns off, an annotation row that drifts once a segment
length changes. Here the construct is defined once as a list of Segments and every diagram
is derived from it, so changing a length re-aligns everything.

Conventions
-----------
* Duplex ladders draw the bottom strand as the plain COMPLEMENT of the top (not the
  reverse complement), written left-to-right and labelled 3'-...-5', so that complementary
  bases share a column. `revcomp()` is provided separately for working out primer sequences.
* Placeholder segments (barcodes, linkers, insert) are NOT complemented -- `AAAAAAAA` is
  barcode A, not eight adenines. This is a property of the Segment, never of the characters.
* Inferred segments render inside <inf>, which the stylesheet gives a dotted underline.
  Colour tags say what a region IS; <inf> says how well it is KNOWN. They nest freely.
"""

from __future__ import annotations

import html
import math
import re
from dataclasses import dataclass, field
from typing import Iterable, Sequence

# --------------------------------------------------------------------------- bases

_COMPLEMENT = str.maketrans("ACGTUacgtuNn", "TGCAAtgcaaNn")

REAL_BASES = set("ACGTUNacgtun")


def complement(seq: str) -> str:
    """Plain complement, preserving order. Non-ACGTN characters pass through unchanged."""
    return seq.translate(_COMPLEMENT)


def revcomp(seq: str) -> str:
    """Reverse complement. Use for primer/landing-site arithmetic, not for drawing ladders."""
    return complement(seq)[::-1]


def is_real_dna(seq: str) -> bool:
    return bool(seq) and all(c in REAL_BASES for c in seq)


def gc_percent(seq: str) -> float:
    s = [c for c in seq.upper() if c in "ACGT"]
    if not s:
        return float("nan")
    return 100.0 * sum(c in "GC" for c in s) / len(s)


# SantaLucia (1998) unified nearest-neighbour parameters: (dH kcal/mol, dS cal/mol/K)
_NN = {
    "AA": (-7.9, -22.2), "TT": (-7.9, -22.2), "AT": (-7.2, -20.4), "TA": (-7.2, -21.3),
    "CA": (-8.5, -22.7), "TG": (-8.5, -22.7), "GT": (-8.4, -22.4), "AC": (-8.4, -22.4),
    "CT": (-7.8, -21.0), "AG": (-7.8, -21.0), "GA": (-8.2, -22.2), "TC": (-8.2, -22.2),
    "CG": (-10.6, -27.2), "GC": (-9.8, -24.4), "GG": (-8.0, -19.9), "CC": (-8.0, -19.9),
}
_R = 1.987


def tm(seq: str, primer_molar: float = 0.5e-6, na_molar: float = 0.05) -> float:
    """Nearest-neighbour melting temperature in degrees C.

    Absolute values depend on buffer (Q5 is not 50 mM Na+); use these comparatively.
    """
    s = seq.upper()
    if not is_real_dna(s) or len(s) < 2:
        raise ValueError(f"tm() needs real DNA of length >= 2, got {seq!r}")
    dh = ds = 0.0
    for i in range(len(s) - 1):
        h, d = _NN[s[i:i + 2]]
        dh += h
        ds += d
    for end in (s[0], s[-1]):                       # helix initiation
        if end in "GC":
            dh += 0.1
            ds += -2.8
        else:
            dh += 2.3
            ds += 4.1
    ds += 0.368 * (len(s) - 1) * math.log(na_molar)  # salt correction
    return (dh * 1000.0) / (ds + _R * math.log(primer_molar / 4.0)) - 273.15


# ------------------------------------------------------------------------ segments

FEATURE_ROLES = (
    "umi",
    "cell_barcode",
    "sample_index",
    "feature_barcode",
    "spatial_barcode",
    "guide_barcode",
    "inline_barcode",
)
FEATURE_ENCODINGS = ("random", "whitelist", "fixed", "combinatorial", "unknown")
_FEATURE_ID = re.compile(r"^[a-z][a-z0-9_.-]*$")


@dataclass(frozen=True)
class MolecularFeature:
    """Machine-readable meaning of an identifier-bearing molecular region.

    This metadata is deliberately independent of ``Segment.tag`` (presentation), the
    displayed bases (which may be placeholders), and any particular interchange format.
    ``id`` is the stable logical identity carried through copying and strand transforms;
    multipart barcodes share ``group`` and use distinct ``part`` values.
    """

    id: str
    role: str
    encoding: str
    group: str = ""
    part: str = ""
    whitelist: str = ""
    note: str = ""

    def __post_init__(self) -> None:
        if not _FEATURE_ID.fullmatch(self.id):
            raise ValueError(
                f"feature id {self.id!r} must start with a lowercase letter and contain "
                "only lowercase letters, digits, '.', '_' or '-'"
            )
        if self.role not in FEATURE_ROLES:
            raise ValueError(f"unknown feature role {self.role!r}; choose from {FEATURE_ROLES}")
        if self.encoding not in FEATURE_ENCODINGS:
            raise ValueError(
                f"unknown feature encoding {self.encoding!r}; choose from {FEATURE_ENCODINGS}"
            )
        if self.encoding == "whitelist" and not self.whitelist:
            raise ValueError(f"whitelist feature {self.id!r} needs a whitelist reference")
        if self.whitelist and self.encoding not in ("whitelist", "combinatorial"):
            raise ValueError(
                f"feature {self.id!r}: a whitelist is incompatible with {self.encoding!r}"
            )
        if self.part and not self.group:
            raise ValueError(f"feature {self.id!r}: multipart feature needs a group")

    @property
    def label(self) -> str:
        return "UMI" if self.role == "umi" else self.role.replace("_", " ")


def feature(identifier: str, role: str, encoding: str, **kw) -> MolecularFeature:
    """Concise constructor used by protocol modules at the point a region is defined."""
    return MolecularFeature(identifier, role, encoding, **kw)


@dataclass
class Segment:
    """One named region of a construct.

    name        label used in annotation rows ('' to leave unlabelled)
    top         top-strand text
    tag         colour element name, e.g. 'p7', 'cbc', 't7'; None for uncoloured
    placeholder stand-in text (barcodes, linkers, insert) -- never complemented
    inferred    wrap in <inf>; the region is a guess, not documented
    bottom      explicit bottom-strand text, overriding the computed complement
    note        free text, surfaced by describe()
    feature     optional format-neutral identifier semantics; preserved through transforms
    """
    name: str
    top: str
    tag: str | None = None
    placeholder: bool = False
    inferred: bool = False
    bottom: str | None = None
    note: str = ""
    feature: MolecularFeature | None = None

    def __post_init__(self) -> None:
        if self.bottom is not None and len(self.bottom) != len(self.top):
            raise ValueError(
                f"segment {self.name!r}: bottom strand length {len(self.bottom)} "
                f"!= top strand length {len(self.top)}"
            )
        if not self.placeholder and self.top and not is_real_dna(self.top):
            raise ValueError(
                f"segment {self.name!r}: {self.top!r} is not real DNA -- "
                f"set placeholder=True if this is a stand-in"
            )

    def __len__(self) -> int:
        return len(self.top)

    def bottom_text(self) -> str:
        if self.bottom is not None:
            return self.bottom
        if self.placeholder:
            if self.is_role_token():
                return self.top
            # lowercase marks "complement of a placeholder": AAAAAAAA (barcode A) -> aaaaaaaa,
            # which can never be mistaken for the dA/dT junction base or for real bases
            return self.top.lower()
        return complement(self.top)

    def is_role_token(self) -> bool:
        """Whether this placeholder is prose for a region, rather than molecular text."""
        return self.placeholder and self.top.startswith("[") and self.top.endswith("]")

    def hover_text(self) -> str | None:
        """Model-derived SVG tooltip for this region.

        A melting temperature is useful for a concrete annealing region, but not for a
        short linker, a placeholder, or an ambiguity-bearing region. The threshold avoids
        presenting physically unhelpful values for two- and three-base structural pieces.
        """
        if not self.name:
            return None
        parts = [self.name]
        seq = self.top.upper()
        if not self.placeholder and len(seq) >= 8 and set(seq) <= set("ACGT"):
            parts.append(f"Tm {tm(seq):.1f} °C")
        return " · ".join(parts)


# ----------------------------------------------------------------------- construct

class Construct:
    """An ordered list of Segments, with offset lookup so nothing is counted by hand."""

    def __init__(self, segments: Sequence[Segment], name: str = "") -> None:
        self.segments = list(segments)
        self.name = name
        seen: set[str] = set()
        features: dict[str, MolecularFeature] = {}
        for s in self.segments:
            if s.name and s.name in seen:
                raise ValueError(f"duplicate segment name {s.name!r} in construct {name!r}")
            if s.name:
                seen.add(s.name)
            if s.feature:
                earlier = features.setdefault(s.feature.id, s.feature)
                if earlier != s.feature:
                    raise ValueError(
                        f"feature id {s.feature.id!r} has conflicting definitions in "
                        f"construct {name!r}"
                    )

    def __len__(self) -> int:
        return sum(len(s) for s in self.segments)

    def __iter__(self):
        return iter(self.segments)

    def top(self) -> str:
        return "".join(s.top for s in self.segments)

    def bottom(self) -> str:
        return "".join(s.bottom_text() for s in self.segments)

    def get(self, name: str) -> Segment:
        for s in self.segments:
            if s.name == name:
                return s
        raise KeyError(f"no segment named {name!r} in construct {self.name!r}")

    def span(self, name: str) -> tuple[int, int]:
        """(start, end) character offsets of a named segment, 0-based half-open."""
        pos = 0
        for s in self.segments:
            if s.name == name:
                return pos, pos + len(s)
            pos += len(s)
        raise KeyError(f"no segment named {name!r} in construct {self.name!r}")

    def offset(self, name: str) -> int:
        return self.span(name)[0]

    def slice(self, first: str | None = None, last: str | None = None,
              name: str = "") -> "Construct":
        """Sub-construct from segment `first` through segment `last`, inclusive."""
        names = [s.name for s in self.segments]
        i = names.index(first) if first else 0
        j = names.index(last) + 1 if last else len(self.segments)
        return Construct(self.segments[i:j], name=name or self.name)

    def describe(self) -> str:
        rows = [f"{'segment':24s} {'start':>6s} {'len':>5s} {'tag':<6s} flags"]
        pos = 0
        for s in self.segments:
            flags = ",".join(
                f for f, on in (("placeholder", s.placeholder), ("inferred", s.inferred)) if on
            )
            if s.feature:
                flags += (("," if flags else "")
                          + f"{s.feature.role}/{s.feature.encoding}:{s.feature.id}")
            rows.append(f"{s.name or '-':24s} {pos:6d} {len(s):5d} {s.tag or '-':<6s} {flags}")
            pos += len(s)
        rows.append(f"{'TOTAL':24s} {'':6s} {pos:5d}")
        return "\n".join(rows)


# ------------------------------------------------------------------------ rendering

def _wrap(text: str, tag: str | None, inferred: bool) -> str:
    out = html.escape(text)
    for t in reversed((tag or "").split("+")):      # "s5+unp" nests both elements
        if t:
            out = f"<{t}>{out}</{t}>"
    if inferred:
        out = f"<inf>{out}</inf>"
    return out


@dataclass
class Row:
    """One rendered line, optionally with geometry derived from molecular semantics.

    ``plain`` and ``html`` retain the character-grid representation for copying and the
    legacy renderer. SVG uses ``visual`` rather than inferring chemistry from punctuation.
    ``chunk_titles`` carries model-derived hover text for corresponding chunks.
    """
    chunks: list[tuple[str, str | None, bool]] = field(default_factory=list)
    indent: int = 0
    prefix: str = ""
    suffix: str = ""
    visual: "StrandVisual | SpanVisual | ArrowVisual | CommentVisual | None" = None
    chunk_titles: list[str | None] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.chunk_titles and len(self.chunk_titles) != len(self.chunks):
            raise ValueError("chunk_titles must be empty or parallel to chunks")

    def plain(self) -> str:
        return " " * self.indent + self.prefix + "".join(c[0] for c in self.chunks) + self.suffix

    def html(self) -> str:
        body = "".join(_wrap(t, g, i) for t, g, i in self.chunks)
        # escape < > & (arrows live here) but leave apostrophes as-is for readability
        esc = lambda t: html.escape(t, quote=False)
        return " " * self.indent + esc(self.prefix) + body + esc(self.suffix)


@dataclass(frozen=True)
class StrandVisual:
    """Arrow-shaped molecular strand spanning character columns ``start:end``."""
    start: int
    end: int
    direction: str

    def __post_init__(self) -> None:
        if self.direction not in ("left", "right") or self.end <= self.start:
            raise ValueError("strand visual needs a non-empty span and left/right direction")


@dataclass(frozen=True)
class SpanVisual:
    """A binding or annealing span located from segment coordinates."""
    start: int
    end: int
    label: str = ""

    def __post_init__(self) -> None:
        if self.end <= self.start:
            raise ValueError("span visual needs a non-empty span")


@dataclass(frozen=True)
class ArrowVisual:
    """Direction of synthesis or enzyme travel, located from a strand end."""
    start: int
    end: int
    direction: str
    label: str = ""

    def __post_init__(self) -> None:
        if self.direction not in ("left", "right") or self.end <= self.start:
            raise ValueError("arrow visual needs a non-empty span and left/right direction")


@dataclass(frozen=True)
class CommentVisual:
    """Prose annotation anchored to a molecular column."""
    start: int
    text: str


@dataclass(frozen=True)
class MolecularState:
    """One rendered molecular state in an ordered reaction workflow."""
    name: str
    rows: tuple[Row, ...]

    def __post_init__(self) -> None:
        if not self.rows:
            raise ValueError(f"molecular state {self.name!r} has no drawing")


@dataclass(frozen=True)
class Reaction:
    """A transition whose input is always the preceding workflow state."""
    action: str
    before: MolecularState
    after: MolecularState
    note: str = ""

    def __post_init__(self) -> None:
        if not self.action.strip():
            raise ValueError("reaction action must not be empty")


class Workflow:
    """Construct a linear reaction path without separately maintained before/after links.

    ``react`` always consumes ``current`` and makes its output the next current state, so
    a renderer cannot accidentally connect a reaction to a stale or unrelated drawing.
    """

    def __init__(self, initial: MolecularState):
        self.initial = initial
        self.reactions: list[Reaction] = []
        self.current = initial

    def react(self, action: str, rows: Iterable[Row], *, name: str = "", note: str = ""):
        after = MolecularState(name, tuple(rows))
        reaction = Reaction(action, self.current, after, note)
        self.reactions.append(reaction)
        self.current = after
        return self

    @property
    def states(self) -> tuple[MolecularState, ...]:
        return (self.initial, *(r.after for r in self.reactions))


def workflow_from_sections(sections, *, initial_rows: Iterable[Row] | None = None,
                           initial_name: str = "Starting material") -> Workflow:
    """Adapt ordered ``(action, resulting rows, note)`` sections to a checked workflow.

    Existing protocol modules can therefore gain reaction rendering centrally. Modules
    with a known molecular input should pass ``initial_rows``; otherwise the renderer is
    explicit that the first input is only identified as starting material.
    """
    start = tuple(initial_rows or (Row(chunks=[("[starting material]", None, False)]),))
    workflow = Workflow(MolecularState(initial_name, start))
    for action, rows, note in sections:
        workflow.react(action, rows, note=note)
    return workflow


def strand_row(con: Construct, strand: str = "top", indent: int = 0,
               prefix: str | None = None, suffix: str | None = None) -> Row:
    """A full duplex strand. Defaults label top as 5'-...-3' and bottom as 3'-...-5'."""
    if strand not in ("top", "bottom"):
        raise ValueError("strand must be 'top' or 'bottom'")
    if prefix is None:
        prefix = "5'- " if strand == "top" else "3'- "
    if suffix is None:
        suffix = " -3'" if strand == "top" else " -5'"
    chunks = [
        ((s.top if strand == "top" else s.bottom_text()), s.tag, s.inferred)
        for s in con if len(s)
    ]
    start = indent + len(prefix)
    return Row(chunks=chunks, indent=indent, prefix=prefix, suffix=suffix,
               visual=StrandVisual(start, start + len(con),
                                   "right" if strand == "top" else "left"),
               chunk_titles=[s.hover_text() for s in con if len(s)])


def annotation_rows(con: Construct, indent: int = 0, prefix_width: int = 4,
                    gap: int = 2) -> list[Row]:
    """Label rows beneath a duplex, each label left-aligned under its segment.

    Labels that would collide are pushed down a level, so this never silently overlaps.
    """
    placed: list[list[tuple[int, str, str | None, bool]]] = []
    pos = 0
    for s in con:
        if s.name:
            for level in placed:
                end = level[-1][0] + len(level[-1][1])
                if pos >= end + gap:
                    level.append((pos, s.name, s.tag, s.inferred))
                    break
            else:
                placed.append([(pos, s.name, s.tag, s.inferred)])
        pos += len(s)

    rows: list[Row] = []
    for level in placed:
        chunks: list[tuple[str, str | None, bool]] = []
        col = 0
        for start, label, tag, inf in level:
            if start > col:
                chunks.append((" " * (start - col), None, False))
            chunks.append((label, tag, inf))
            col = start + len(label)
        rows.append(Row(chunks=chunks, indent=indent + prefix_width))
    return rows


def feature_rows(con: Construct, indent: int = 0, prefix_width: int = 0) -> list[Row]:
    """Rows locating identifier semantics beneath a construct.

    One row per occurrence keeps multipart barcodes readable and makes the marker's
    starting column exact.  It is intentionally derived from ``Segment.feature`` rather
    than colour tags or prose names.
    """
    rows, pos = [], 0
    for s in con:
        if s.feature:
            f = s.feature
            details = [f.label]
            if f.part:
                details.append(f.part)
            rows.append(Row(chunks=[(" " * pos + "^ " + " · ".join(details),
                                     s.tag, s.inferred)],
                            indent=indent + prefix_width))
        pos += len(s)
    return rows


def duplex_rows(con: Construct, label: str = "", *,
                bottom: Sequence[Segment] | None = None,
                unpaired: Sequence[str] = ()) -> list[Row]:
    """A complete final-library duplex with aligned names and feature markers.

    ``Scene`` owns the molecular placement, while the two annotation helpers operate in
    construct coordinates.  Keeping the gutter calculation here prevents callers from
    guessing how much room the scene's strand label and 5'/3' end label consume.
    """
    scene = Scene.duplex(list(con), label=label, bottom=bottom, unpaired=unpaired)
    sequence_column = (len(label) + 1 if label else 0) + len("5'- ")
    return [*scene.rows(),
            *annotation_rows(con, prefix_width=sequence_column),
            *feature_rows(con, prefix_width=sequence_column)]


def junction_row(con: Construct, left: str, right: str, text: str = "ligation",
                 ch: str = "*", indent: int = 0, prefix_width: int = 4) -> Row:
    """Mark the exact boundary between two adjacent construct segments.

    Naming both sides makes the marker follow the chemistry when segment lengths change,
    and refuses to draw it if another segment is inserted at the claimed junction.
    """
    names = [s.name for s in con]
    try:
        i, j = names.index(left), names.index(right)
    except ValueError as exc:
        raise ValueError(f"ligation junction needs segments {left!r} and {right!r}") from exc
    if j != i + 1:
        raise ValueError(f"ligation junction segments are not adjacent: {left!r}, {right!r}")
    boundary = con.span(left)[1]
    marker = ch * 2 + (" " + text if text else "")
    return Row(chunks=[(" " * (boundary - 1) + marker, None, False)],
               indent=indent + prefix_width)


def primer_row(con: Construct, seq_segments: Sequence[Segment], anchor: str,
               direction: str = ">", indent: int = 0, prefix_width: int = 4,
               arrow: int = 8, label_5p: bool = True) -> Row:
    """A primer drawn above/below the construct, positioned by segment name.

    `anchor` is the segment the primer's 5'-most drawn base sits over (for '>'), so the
    column arithmetic is done here rather than by counting spaces in a heredoc.
    """
    start = con.offset(anchor)
    chunks: list[tuple[str, str | None, bool]] = []
    pre = "5'- " if label_5p else ""
    pad = start + prefix_width - len(pre)
    if pad < 0:
        raise ValueError(
            f"primer anchored at {anchor!r} starts before column 0 "
            f"(needs {-pad} more columns of indent)"
        )
    chunks.append((" " * pad, None, False))
    if pre:
        chunks.append((pre, None, False))
    for s in seq_segments:
        chunks.append((s.top, s.tag, s.inferred))
    chunks.append(("-" * arrow + (">" if direction == ">" else ""), None, False))
    return Row(chunks=chunks, indent=indent)


def panel_legacy(rows: Iterable[Row], cls: str = "long",
                 caption: str | None = None) -> str:
    """The original character-grid ``<pre>`` renderer, retained as a fallback."""
    body = []
    if caption:
        body.append(f"<i>{html.escape(caption)}</i>")
        body.append("")
    body.extend(r.html() for r in rows)
    return f'<pre>\n<align class="{cls}">\n' + "\n".join(body) + "\n</align>\n</pre>"


def _svg_classes(tag: str | None, inferred: bool) -> str:
    classes = [f"chem-{t}" for t in (tag or "").split("+") if t]
    if inferred:
        classes.append("chem-inferred")
    return " ".join(classes)


def _svg_row(row: Row, *, x: float, y: float, font: float, cell: float) -> list[str]:
    """Render one row; molecular geometry comes from typed metadata, never ASCII art."""
    baseline = y + font
    visual = row.visual
    out: list[str] = []

    if isinstance(visual, StrandVisual):
        x1, x2 = x + visual.start * cell - 4, x + visual.end * cell + 4
        top, bottom, middle = y - 1, y + font + 4, y + (font + 3) / 2
        head = min(10.0, max(5.0, (x2 - x1) / 4))
        if visual.direction == "right":
            points = ((x1, top), (x2 - head, top), (x2, middle),
                      (x2 - head, bottom), (x1, bottom))
        else:
            points = ((x2, top), (x1 + head, top), (x1, middle),
                      (x1 + head, bottom), (x2, bottom))
        coords = " ".join(f"{a:.2f},{b:.2f}" for a, b in points)
        out.append(f'<polygon class="chem-strand-box" points="{coords}"/>')

    if isinstance(visual, SpanVisual):
        x1, x2 = x + visual.start * cell, x + visual.end * cell
        rule_y = y + 3.0
        out.append(f'<path class="chem-binding-span" d="M {x1:.2f} {rule_y + 4:.2f} '
                   f'V {rule_y:.2f} H {x2:.2f} V {rule_y + 4:.2f}"/>')
        if visual.label:
            out.append(f'<text class="chem-binding-label" x="{x2 + 7:.2f}" '
                       f'y="{baseline:.2f}">{html.escape(visual.label)}</text>')
        return out

    if isinstance(visual, ArrowVisual):
        x1, x2 = x + visual.start * cell, x + visual.end * cell
        line_y = y + font * .58
        tip = x1 if visual.direction == "left" else x2
        tail = x2 if visual.direction == "left" else x1
        out.append(f'<line class="chem-process-arrow" x1="{tail:.2f}" y1="{line_y:.2f}" '
                   f'x2="{tip:.2f}" y2="{line_y:.2f}"/>')
        sign = 1 if visual.direction == "left" else -1
        points = ((tip, line_y), (tip + sign * 7, line_y - 4),
                  (tip + sign * 7, line_y + 4))
        coords = " ".join(f"{a:.2f},{b:.2f}" for a, b in points)
        out.append(f'<polygon class="chem-process-arrowhead" points="{coords}"/>')
        if visual.label:
            label_x = x1 - 7 if visual.direction == "left" else x2 + 7
            anchor = "end" if visual.direction == "left" else "start"
            out.append(f'<text class="chem-process-label" x="{label_x:.2f}" '
                       f'y="{baseline:.2f}" text-anchor="{anchor}">'
                       f'{html.escape(visual.label)}</text>')
        return out

    if isinstance(visual, CommentVisual):
        out.append(f'<text class="chem-comment" x="{x + visual.start * cell:.2f}" '
                   f'y="{baseline:.2f}">{html.escape(visual.text)}</text>')
        return out

    start = " " * row.indent + row.prefix
    pieces = [html.escape(start, quote=False)]
    titles = row.chunk_titles or [None] * len(row.chunks)
    for (text, tag, inferred), title in zip(row.chunks, titles):
        classes = _svg_classes(tag, inferred)
        if title:
            classes = (classes + " chem-has-tip").strip()
        attr = f' class="{classes}"' if classes else ""
        tooltip = f"<title>{html.escape(title)}</title>" if title else ""
        pieces.append(f"<tspan{attr}>{tooltip}{html.escape(text, quote=False)}</tspan>")
    pieces.append(html.escape(row.suffix, quote=False))
    out.append(f'<text x="{x:g}" y="{baseline:.2f}" xml:space="preserve">'
               + "".join(pieces) + "</text>")
    return out


def panel_svg(rows: Iterable[Row], cls: str = "long",
              caption: str | None = None) -> str:
    """Render a selectable, fixed-scale SVG character grid.

    Every row remains one SVG text node, with coloured sequence regions as ``tspan``
    children.  Decorative geometry can therefore evolve without turning bases into
    paths or compromising copy/paste.  The intrinsic pixel width deliberately follows
    the longest row; CSS scrolls it instead of shrinking the font.
    """
    rows = list(rows)
    font = 11.8 if cls == "long" else 13.4
    cell = font * 0.602                 # IBM Plex Mono advance width
    line = font * 1.42
    pad_x, pad_y = 14.0, 12.0
    columns = max([len(r.plain()) for r in rows] or [1])
    width = max(280.0, pad_x * 2 + columns * cell)
    height = max(42.0, pad_y * 2 + max(1, len(rows)) * line)
    text_rows = []
    for i, row in enumerate(rows):
        text_rows.extend(_svg_row(row, x=pad_x, y=pad_y + i * line,
                                  font=font, cell=cell))
    label = html.escape(caption or "Molecular construct diagram", quote=True)
    cap = (f'<figcaption>{html.escape(caption)}</figcaption>' if caption else "")
    svg = (f'<svg class="chem-svg {html.escape(cls)}" width="{width:.0f}" '
           f'height="{height:.0f}" viewBox="0 0 {width:.2f} {height:.2f}" '
           f'role="img" aria-label="{label}">\n' + "\n".join(text_rows) + "\n</svg>")
    return f'<figure class="chem-panel">{cap}<div class="diagram-scroll">{svg}</div></figure>'


def _svg_rows(rows: Sequence[Row], *, x: float, y: float, font: float,
              line: float) -> tuple[list[str], float]:
    """SVG elements for one state, returning elements and the next y coordinate."""
    out = []
    cell = font * 0.602
    for row in rows:
        out.extend(_svg_row(row, x=x, y=y, font=font, cell=cell))
        y += line
    return out, y


def workflow_panel(workflow: Workflow, cls: str = "long",
                   caption: str = "Reaction workflow") -> str:
    """Render molecular states once, connected by real vector reaction arrows.

    The SVG is derived only from ``Workflow``. Protocol pages do not draw arrows or copy
    input states themselves, and the workflow constructor guarantees every arrow starts
    at the immediately preceding molecular state.
    """
    font = 11.8 if cls == "long" else 13.4
    cell = font * 0.602
    line = font * 1.42
    pad_x, pad_y = 14.0, 12.0
    arrow_h, label_h = 54.0, 22.0
    states = workflow.states
    row_columns = max((len(row.plain()) for state in states for row in state.rows), default=1)
    action_columns = max((len(r.action) for r in workflow.reactions), default=1) + 13
    width = max(420.0, pad_x * 2 + max(row_columns, action_columns) * cell)
    x_arrow = pad_x + 13.0
    y = pad_y
    elements = []

    if workflow.initial.name:
        elements.append(f'<text class="chem-state-label" x="{pad_x:g}" y="{y + font:.2f}">'
                        f'{html.escape(workflow.initial.name)}</text>')
        y += line
    block, y = _svg_rows(workflow.initial.rows, x=pad_x, y=y, font=font, line=line)
    elements.extend(block)

    for reaction in workflow.reactions:
        y += 7.0
        top, bottom = y, y + arrow_h
        elements.append(f'<line class="chem-reaction-arrow" x1="{x_arrow:g}" y1="{top:.2f}" '
                        f'x2="{x_arrow:g}" y2="{bottom - 8:.2f}"/>')
        elements.append(f'<path class="chem-reaction-arrowhead" d="M {x_arrow - 5:g} '
                        f'{bottom - 10:.2f} L {x_arrow:g} {bottom:.2f} L {x_arrow + 5:g} '
                        f'{bottom - 10:.2f} Z"/>')
        label_width = max(70.0, len(reaction.action) * cell + 16.0)
        label_y = top + (arrow_h - label_h) / 2
        elements.append(f'<rect class="chem-reaction-box" x="{x_arrow + 14:g}" '
                        f'y="{label_y:.2f}" width="{label_width:.2f}" height="{label_h:g}" '
                        f'rx="6"/>')
        elements.append(f'<text class="chem-reaction-label" x="{x_arrow + 22:g}" '
                        f'y="{label_y + font + 3:.2f}">{html.escape(reaction.action)}</text>')
        y = bottom + 5.0
        block, y = _svg_rows(reaction.after.rows, x=pad_x, y=y, font=font, line=line)
        elements.extend(block)

    height = max(42.0, y + pad_y)
    label = html.escape(caption, quote=True)
    notes = [r for r in workflow.reactions if r.note]
    note_html = ""
    if notes:
        items = "".join(f'<li><b>{html.escape(r.action)}.</b> {html.escape(r.note)}</li>'
                        for r in notes)
        note_html = f'<figcaption><ol class="reaction-notes">{items}</ol></figcaption>'
    svg = (f'<svg class="chem-svg chem-workflow {html.escape(cls)}" width="{width:.0f}" '
           f'height="{height:.0f}" viewBox="0 0 {width:.2f} {height:.2f}" role="img" '
           f'aria-label="{label}">\n' + "\n".join(elements) + "\n</svg>")
    return (f'<figure class="chem-panel reaction-panel"><div class="diagram-scroll">'
            f'{svg}</div>{note_html}</figure>')


def panel(rows: Iterable[Row], cls: str = "long", caption: str | None = None,
          renderer: str = "svg") -> str:
    """Render a diagram; SVG is default, ``renderer="legacy"`` keeps the old form."""
    rows = list(rows)
    if renderer == "svg":
        return panel_svg(rows, cls, caption)
    if renderer == "legacy":
        return panel_legacy(rows, cls, caption)
    raise ValueError("panel renderer must be 'svg' or 'legacy'")


def oligo(name: str, segments: Sequence[Segment], five: str = "5'-", three: str = "-3'",
          mods: str = "") -> str:
    """One raw-text line for ordering or copying an oligo sequence."""
    body = "".join(_wrap(s.top, s.tag, s.inferred) for s in segments)
    lead = f"{five} {mods} " if mods else f"{five} "
    return (f"<p>{html.escape(name)}: {html.escape(lead, quote=False)}"
            f"{body} {html.escape(three, quote=False)}</p>")


# ------------------------------------------------------------------- loops
def bridge(left: int, right: int, inner: Sequence[Sequence[tuple]] = (),
           label: str = "", riser: int = 1) -> list[Row]:
    """Draw a loop joining two columns, with content carried inside it.

    Used where a molecule leaves the duplex, travels, and comes back -- a padlock
    probe's backbone arching from one annealed arm to the other, a hairpin, a lariat.
    `left` and `right` are the columns the two risers come down on, so the caller can
    anchor them to real sequence positions rather than eyeballing the width.

    Returns rows top-down: the arc, then one row per `inner` line, then `riser` plain
    rows. Each `inner` line is a sequence of (text, tag, inferred) chunks.
    """
    width = right - left + 1
    if width < 4:
        raise ValueError(f"bridge needs at least 4 columns, got {width}")
    body = width - 2
    if label:
        lab = f" {label} "
        if len(lab) > body:
            lab = lab[:body]
        pad = body - len(lab)
        bar = "-" * (pad // 2) + lab + "-" * (pad - pad // 2)
    else:
        bar = "-" * body
    rows = [Row(chunks=[(" " * left + "." + bar + ".", None, False)])]
    for line in inner:
        used = sum(len(c[0]) for c in line)
        if used > body:
            raise ValueError(f"bridge content {used} wide does not fit in {body}")
        lead = (body - used) // 2
        rows.append(Row(chunks=[(" " * left + "|" + " " * lead, None, False), *line,
                                (" " * (body - used - lead) + "|", None, False)]))
    for _ in range(riser):
        rows.append(Row(chunks=[(" " * left + "|" + " " * body + "|", None, False)]))
    return rows


def circle_rows(con: Construct, closure_label: str = "covalently closed") -> list[Row]:
    """Draw one circular strand opened at its last-to-first bond.

    The sequence stays linear and readable, while the return path makes the topology
    explicit.  Because the displayed break is always the boundary between the final and
    first segments, callers cannot accidentally mark some unrelated internal bond as the
    circularisation junction.
    """
    if not len(con):
        raise ValueError("a circle needs at least one base")
    lead = "  .-> "
    molecule = strand_row(con, prefix=lead, suffix=" -.")
    left, right = lead.index("."), len(molecule.plain()) - 1
    width = right - left - 1
    label = (f" ** {closure_label}: {con.segments[-1].name} -> "
             f"{con.segments[0].name} ** ")
    if len(label) > width:
        label = " ** circularisation ** "
    pad = width - len(label)
    close = (" " * left + "'" + "-" * (pad // 2) + label
             + "-" * (pad - pad // 2) + "'")
    return [molecule, Row(chunks=[(close, None, False)])]


# ------------------------------------------------------------------- pairing
# One rule for "may these two characters be drawn in the same column of a duplex", used both
# by Scene (refuses to build a wrong drawing) and by tools/lint_pairing.py (finds wrong
# drawings in already-generated pages).

_IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "U": "T", "N": "ACGT", "W": "AT",
          "S": "CG", "R": "AG", "Y": "CT", "K": "GT", "M": "AC", "B": "CGT", "D": "AGT",
          "H": "ACT", "V": "ACG"}
_WC = {"A": "T", "T": "A", "C": "G", "G": "C"}


def bases_pair(a: str, b: str) -> bool:
    """True if `a` may sit opposite `b` in a drawn duplex.

    Real bases (and IUPAC codes) must be Watson-Crick compatible. Placeholders follow the
    repo convention: the complement of placeholder `X` is drawn `x`, so the same letter in
    opposite case pairs. `.` pairs with `.` (the `XXX...XXX` ellipsis).
    """
    if a == "." or b == ".":
        return a == b
    if a in "Xx" and b in "XxNn":
        return True                       # unknown DNA; some pages draw it X on both strands
    if b in "Xx" and a in "Nn":
        return True
    if a.isalpha() and b.isalpha() and a != b and a.upper() == b.upper():
        return True
    sa, sb = _IUPAC.get(a.upper()), _IUPAC.get(b.upper())
    if sa is None or sb is None:
        return False
    return any(_WC[x] in sb for x in sa)


# --------------------------------------------------------------------- scenes
# Hand-written indents are how drawings go wrong: an oligo-dT two columns into the cDNA, a
# GGG nowhere near its CCC. A Scene takes every strand 5'->3' as written on the order sheet,
# places it by saying WHICH SEGMENT PAIRS WITH WHICH, computes all columns itself, and
# raises if any drawn column does not base-pair. There is no indent argument anywhere.

@dataclass
class _Strand:
    name: str
    segs: list[Segment]          # 5'->3', as ordered
    col: int                     # column of the leftmost DRAWN base
    rev: bool                    # True: drawn 3'->5' left to right
    label: str
    unpaired: tuple[str, ...] = ()
    partner: str | None = None
    mod5: str = ""               # e.g. "p" -> drawn 5'-p / p-5'
    mod3: str = ""

    def ends(self) -> tuple[str, str]:
        """(left label, right label) as drawn."""
        f5 = "5'-" + self.mod5 + " " if not self.rev else " " + self.mod5 + "-5'"
        f3 = " -3'" + (" " + self.mod3 if self.mod3 else "") if not self.rev else \
            (self.mod3 + " " if self.mod3 else "") + "3'- "
        return (f5, f3) if not self.rev else (f3, f5)

    def drawn(self) -> list[Segment]:
        if not self.rev:
            return self.segs
        return [Segment(s.name, s.top if s.is_role_token() else s.top[::-1],
                        s.tag, s.placeholder, s.inferred,
                        None if s.bottom is None else
                        (s.bottom if s.is_role_token() else s.bottom[::-1]), s.note,
                        s.feature)
                for s in reversed(self.segs)]

    def text(self) -> str:
        return "".join(s.top for s in self.drawn())

    def span(self, seg: str) -> tuple[int, int]:
        pos = self.col
        for s in self.drawn():
            if s.name == seg:
                return pos, pos + len(s)
            pos += len(s)
        raise KeyError(f"strand {self.name!r} has no segment {seg!r}")

    def end(self) -> int:
        return self.col + len(self.text())

    def cells(self) -> dict[int, tuple[str, str]]:
        """column -> (char, segment name)"""
        out, pos = {}, self.col
        for s in self.drawn():
            for ch in s.top:
                out[pos] = (ch, s.name)
                pos += 1
        return out


def complement_segments(segs: Sequence[Segment], suffix: str = "'") -> list[Segment]:
    """The strand that pairs with `segs` (given 5'->3'), itself returned 5'->3'.

    Real bases are reverse-complemented; placeholders become the reversed lowercase stand-in
    (X -> x), the same convention the duplex ladders use. A segment with an explicit
    `bottom` (e.g. a T:A junction placeholder) contributes exactly that bottom text.
    """
    out = []
    for s in reversed(segs):
        top = s.bottom_text() if s.is_role_token() else s.bottom_text()[::-1]
        out.append(Segment(s.name + suffix if s.name else "", top, s.tag,
                           s.placeholder, s.inferred, feature=s.feature))
    return out


class Scene:
    """Strands placed by pairing, never by indent. Render with `rows()`.

        sc = Scene()
        sc.strand("mRNA", mrna_segs)                                # first strand: column 0
        sc.anneal("oligo-dT", dt_segs, to="mRNA", pair=("dT", "polyA"), shift=6)
        sc.mark("mRNA", "polyA", "anywhere in the tract")

    `anneal` draws the new strand antiparallel to its partner, with segment `pair[0]`'s
    drawn left edge `shift` columns right of `pair[1]`'s, then checks EVERY overlapping
    column with `bases_pair`. Segments listed in `unpaired` (e.g. a mismatched 5' tail)
    are exempt. A drawing that would put a T over an X cannot be produced.
    """

    def __init__(self) -> None:
        self.strands: dict[str, _Strand] = {}
        self.order: list[tuple[str, str]] = []      # ("strand"|"mark"|"arrow"|"blank", key)
        self.extras: dict[str, tuple] = {}
        self._same_line: list[tuple[str, str]] = []
        self._footers: set[str] = set()

    @classmethod
    def duplex(cls, top: Sequence[Segment], on: str | None = None, label: str = "",
               bottom: Sequence[Segment] | None = None, **kw) -> "Scene":
        """Strand "top" plus its pairing strand "bottom" (default: the full complement),
        paired on segment `on` (default: the first). kw (mod5=, unpaired=...) go to anneal."""
        on = top[0].name if on is None else on
        bottom = complement_segments(top) if bottom is None else bottom
        sc = cls()
        sc.strand("top", top, label=label)
        sc.anneal("bottom", bottom, to="top", pair=(on + "'", on), label=label, **kw)
        return sc

    # -- placement
    def strand(self, name: str, segs: Sequence[Segment], label: str | None = None,
               rev: bool = False, at: int = 0, mod5: str = "", mod3: str = "") -> str:
        """A free strand (no partner). Normally only the first one in a scene."""
        self._add(_Strand(name, list(segs), at, rev, name if label is None else label,
                          mod5=mod5, mod3=mod3))
        return name

    def anneal(self, name: str, segs: Sequence[Segment], to: str, pair: tuple[str, str],
               shift: int = 0, label: str | None = None, above: bool = False,
               unpaired: Sequence[str] = (), mod5: str = "", mod3: str = "") -> str:
        p = self.strands[to]
        st = _Strand(name, list(segs), 0, not p.rev, name if label is None else label,
                     tuple(unpaired), to, mod5, mod3)
        mine, theirs = pair
        st.col = p.span(theirs)[0] + shift - (st.span(mine)[0] - st.col)
        self._verify(st, p)
        self._add(st, before=to if above else None, after=None if above else to)
        return name

    def _verify(self, a: _Strand, b: _Strand) -> None:
        ca, cb = a.cells(), b.cells()
        both = sorted(set(ca) & set(cb))
        if not both:
            raise ValueError(f"{a.name!r} is drawn against {b.name!r} but no column overlaps")
        bad = [c for c in both
               if ca[c][1] not in a.unpaired and cb[c][1] not in b.unpaired
               and not bases_pair(ca[c][0], cb[c][0])]
        if bad:
            c = bad[0]
            raise ValueError(
                f"{a.name!r} vs {b.name!r}: {len(bad)} drawn column(s) do not pair, first at "
                f"{a.name}:{ca[c][1]} {ca[c][0]!r} over {b.name}:{cb[c][1]} {cb[c][0]!r}\n"
                f"  {b.text()!r} @ {b.col}\n  {a.text()!r} @ {a.col}")

    def _add(self, st: _Strand, before: str | None = None, after: str | None = None) -> None:
        if st.name in self.strands:
            raise ValueError(f"duplicate strand {st.name!r}")
        self.strands[st.name] = st
        item = ("strand", st.name)
        keys = [k for _, k in self.order]
        if before is not None:
            self.order.insert(keys.index(before), item)
        elif after is not None:
            i = keys.index(after) + 1
            while i < len(self.order) and self.order[i][0] in ("mark", "arrow"):
                i += 1                                  # keep a strand's marks attached to it
            self.order.insert(i, item)
        else:
            self.order.append(item)

    # -- decorations, all positioned from segments
    def mark(self, strand: str, seg: str, text: str, ch: str = "^",
             through: str | None = None) -> None:
        """Underline segment `seg` (or `seg`..`through`, in drawn order) and label it."""
        st = self.strands[strand]
        a = st.span(seg)
        b = st.span(through) if through else a
        s, e = min(a[0], b[0]), max(a[1], b[1])
        self._decor("mark", strand,
                    (s, ch * (e - s) + (" " + text if text else ""),
                     SpanVisual(s, e, text)))

    def junction(self, strand: str, left: str, right: str, text: str = "ligation",
                 ch: str = "*") -> None:
        """Mark the exact boundary between two adjacent segments on a drawn strand."""
        st = self.strands[strand]
        drawn = [s.name for s in st.drawn()]
        try:
            i, j = drawn.index(left), drawn.index(right)
        except ValueError as exc:
            raise ValueError(
                f"ligation junction needs {left!r} and {right!r} on strand {strand!r}"
            ) from exc
        if abs(i - j) != 1:
            raise ValueError(
                f"ligation junction segments are not adjacent on {strand!r}: "
                f"{left!r}, {right!r}"
            )
        a, b = sorted((st.span(left), st.span(right)))
        if a[1] != b[0]:
            raise ValueError(f"ligation junction has a gap on strand {strand!r}")
        marker = ch * 2 + (" " + text if text else "")
        self._decor("mark", strand, (a[1] - 1, marker, None))

    def note(self, strand: str, text: str) -> None:
        """A free line under the strand, starting at its first drawn base."""
        start = self.strands[strand].col
        self._decor("mark", strand, (start, text, CommentVisual(start, text)))

    def footer(self, text: str, strand: str, seg: str | None = None) -> None:
        """A line at the very bottom of the scene, starting at segment `seg` of `strand`
        (or the strand's first base). Unlike mark/note it never splits a duplex."""
        st = self.strands[strand]
        col = st.span(seg)[0] if seg else st.col
        key = f"mark{len(self.extras)}"
        self.extras[key] = (col, text, None)
        self.order.append(("mark", key))
        self._footers.add(key)

    def labels(self, strand: str, gap: int = 2) -> None:
        """Segment names as footer lines, each starting on its segment in `strand`, packed
        onto as few lines as possible (the Scene analogue of annotation_rows)."""
        st = self.strands[strand]
        lines: list[str] = []
        for s in st.drawn():
            if not s.name or not s.top:
                continue
            c = st.span(s.name)[0] - st.col
            for i, ln in enumerate(lines):
                if c >= len(ln) + gap:
                    lines[i] = ln.ljust(c) + s.name
                    break
            else:
                lines.append(" " * c + s.name)
        for ln in lines:
            self.footer(ln, strand)

    def stack(self, *names: str) -> None:
        """Reorder strands top-down as `names` (every strand exactly once); each strand keeps
        its own marks and arrows under it. Use when call order cannot express the layout,
        e.g. one primer above the top strand and another below the bottom strand."""
        groups: dict[str, list] = {}
        head: list = []
        cur = None
        feet = [it for it in self.order if it[1] in self._footers]
        for item in self.order:
            if item[1] in self._footers:
                continue
            if item[0] == "strand":
                cur = item[1]
                groups[cur] = []
            (groups[cur] if cur else head).append(item)
        if sorted(names) != sorted(groups):
            raise ValueError(f"stack() needs every strand exactly once: {sorted(groups)}")
        self.order = head + [it for n in names for it in groups[n]] + feet

    def same_line(self, host: str, guest: str) -> None:
        """Draw two collinear but unjoined strands on one line (e.g. a Tn5 non-transferred
        strand 9 nt beyond a fragment's 3' end). Refuses at render if they would collide."""
        self._same_line.append((host, guest))

    def blank(self, before: str) -> None:
        """An empty line directly above strand `before`."""
        keys = [k for _, k in self.order]
        self.order.insert(keys.index(before), ("blank", f"blank{len(self.order)}"))

    def origin(self) -> tuple[int, int]:
        """(screen column of scene column 0, gutter width), as rows() lays it out -- for
        helpers that add rows aligned to the scene (e.g. a padlock bridge)."""
        lefts = [st.col - len(st.ends()[0]) for st in self.strands.values()]
        lefts += [value[0] for value in self.extras.values()]
        gutter = max(len(st.label) for st in self.strands.values()) + 1
        if not any(st.label for st in self.strands.values()):
            gutter = 0
        return gutter - min(lefts), gutter

    def arrow(self, strand: str, text: str, length: int = 8) -> None:
        """Extension arrow off the strand's 3' end, pointing the way synthesis goes."""
        st = self.strands[strand]
        if st.rev:
            body = f"{text} <" + "-" * length
            self._decor("arrow", strand,
                        (st.col - len(body), body,
                         ArrowVisual(st.col - length, st.col, "left", text)))
        else:
            self._decor("arrow", strand,
                        (st.end(), "-" * length + "> " + text,
                         ArrowVisual(st.end(), st.end() + length, "right", text)))

    def _decor(self, kind: str, strand: str, payload: tuple) -> None:
        key = f"{kind}{len(self.extras)}"
        self.extras[key] = payload
        # under the strand, but never between it and a strand annealed directly below it
        keys = [k for _, k in self.order]
        i = keys.index(strand) + 1
        while i < len(self.order):
            k, n = self.order[i]
            if k in ("mark", "arrow") or (k == "strand" and self.strands[n].partner == strand):
                i += 1
            else:
                break
        self.order.insert(i, (kind, key))

    # -- rendering
    def rows(self, omit: Sequence[str] = ()) -> list[Row]:
        """Rendered rows. Strands named in `omit` are left out (their marks stay) -- for
        helpers that redraw those strands themselves, e.g. two padlock arms on one row."""
        shift, gutter = self.origin()
        out: list[Row] = []
        keys: list[str | None] = []
        for kind, key in self.order:
            if kind == "strand":
                st = self.strands[key]
                l_end, r_end = st.ends()
                lead = st.label.ljust(gutter) + " " * (shift - gutter + st.col - len(l_end))
                drawn = [s for s in st.drawn() if s.top]
                seq_start = len(lead) + len(l_end)
                seq_len = sum(len(s.top) for s in drawn)
                out.append(Row(chunks=[(lead + l_end, None, False)]
                               + [(s.top, (s.tag or "") + ("+unp" if s.name in st.unpaired
                                                            else ""), s.inferred)
                                  for s in drawn]
                               + [(r_end, None, False)],
                               visual=StrandVisual(seq_start, seq_start + seq_len,
                                                   "left" if st.rev else "right"),
                               chunk_titles=[None]
                               + [s.hover_text() for s in drawn] + [None]))
                keys.append(key)
            elif kind == "blank":
                out.append(Row())
                keys.append(None)
            else:
                c, text, visual = self.extras[key]
                if isinstance(visual, SpanVisual):
                    visual = SpanVisual(shift + visual.start, shift + visual.end, visual.label)
                elif isinstance(visual, ArrowVisual):
                    visual = ArrowVisual(shift + visual.start, shift + visual.end,
                                         visual.direction, visual.label)
                elif isinstance(visual, CommentVisual):
                    visual = CommentVisual(shift + visual.start, visual.text)
                out.append(Row(chunks=[(" " * (shift + c) + text, None, False)],
                               visual=visual))
                keys.append(None)
        for host, guest in self._same_line:
            hi, gi = keys.index(host), keys.index(guest)
            if self.strands[host].col > self.strands[guest].col:
                hi, gi = gi, hi
            left = out[hi].plain().rstrip()
            lead = out[gi].chunks[0][0]
            if lead[:len(left) + 1].strip():
                raise ValueError(f"{host!r} and {guest!r} overlap on one line")
            h = list(out[hi].chunks)
            h[-1] = (h[-1][0].rstrip(), h[-1][1], h[-1][2])
            out[hi] = Row(chunks=h + [(lead[len(left):], None, False)] + list(out[gi].chunks[1:]))
            out[gi] = None
        drop = {i for i, k in enumerate(keys) if k is not None and k in omit}
        return [r for i, r in enumerate(out) if r is not None and i not in drop]
