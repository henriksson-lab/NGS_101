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
from dataclasses import dataclass, field
from typing import Iterable, Sequence

# --------------------------------------------------------------------------- bases

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

REAL_BASES = set("ACGTNacgtn")


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
    """
    name: str
    top: str
    tag: str | None = None
    placeholder: bool = False
    inferred: bool = False
    bottom: str | None = None
    note: str = ""

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
            # lowercase marks "complement of a placeholder": AAAAAAAA (barcode A) -> aaaaaaaa,
            # which can never be mistaken for the dA/dT junction base or for real bases
            return self.top.lower()
        return complement(self.top)


# ----------------------------------------------------------------------- construct

class Construct:
    """An ordered list of Segments, with offset lookup so nothing is counted by hand."""

    def __init__(self, segments: Sequence[Segment], name: str = "") -> None:
        self.segments = list(segments)
        self.name = name
        seen: set[str] = set()
        for s in self.segments:
            if s.name and s.name in seen:
                raise ValueError(f"duplicate segment name {s.name!r} in construct {name!r}")
            if s.name:
                seen.add(s.name)

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
            rows.append(f"{s.name or '-':24s} {pos:6d} {len(s):5d} {s.tag or '-':<6s} {flags}")
            pos += len(s)
        rows.append(f"{'TOTAL':24s} {'':6s} {pos:5d}")
        return "\n".join(rows)


# ------------------------------------------------------------------------ rendering

def _wrap(text: str, tag: str | None, inferred: bool) -> str:
    out = html.escape(text)
    if tag:
        out = f"<{tag}>{out}</{tag}>"
    if inferred:
        out = f"<inf>{out}</inf>"
    return out


@dataclass
class Row:
    """One rendered line. `indent` is in characters; chunks are (text, tag, inferred)."""
    chunks: list[tuple[str, str | None, bool]] = field(default_factory=list)
    indent: int = 0
    prefix: str = ""
    suffix: str = ""

    def plain(self) -> str:
        return " " * self.indent + self.prefix + "".join(c[0] for c in self.chunks) + self.suffix

    def html(self) -> str:
        body = "".join(_wrap(t, g, i) for t, g, i in self.chunks)
        # escape < > & (arrows live here) but leave apostrophes as-is for readability
        esc = lambda t: html.escape(t, quote=False)
        return " " * self.indent + esc(self.prefix) + body + esc(self.suffix)


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
    return Row(chunks=chunks, indent=indent, prefix=prefix, suffix=suffix)


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


def panel(rows: Iterable[Row], cls: str = "long", caption: str | None = None) -> str:
    """Render rows as one <pre><align> block, with an optional <i> caption inside."""
    body = []
    if caption:
        body.append(f"<i>{html.escape(caption)}</i>")
        body.append("")
    body.extend(r.html() for r in rows)
    return f'<pre>\n<align class="{cls}">\n' + "\n".join(body) + "\n</align>\n</pre>"


def oligo(name: str, segments: Sequence[Segment], five: str = "5'-", three: str = "-3'",
          mods: str = "") -> str:
    """One line of the 'Adapter and primer sequences' list."""
    body = "".join(_wrap(s.top, s.tag, s.inferred) for s in segments)
    lead = f"{five} {mods} " if mods else f"{five} "
    return (f"<p>{html.escape(name)}: {html.escape(lead, quote=False)}"
            f"{body} {html.escape(three, quote=False)}</p>")
