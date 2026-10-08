"""
Sequencing primers (Read 1, Index 1, Index 2, Read 2), and where they land on a library.

One source of truth. The primer SEQUENCES live in illumina.py / nextera.py and nowhere
else; this module only names the standard sets and references those constants. A protocol
declares which primers it uses -- by reference, never by pasting sequence -- and
`section()` renders them for its page with everything else COMPUTED from the protocol's
final library: which strand each primer anneals to, which segments it covers, and the
first bases each read reports. A declared primer that has no exact site raises, unless the
mismatch is declared with `expect_mismatch="why"`.
"""

from __future__ import annotations

import html
from dataclasses import dataclass, replace

import illumina as il
import nextera as nx
from chemdraw import Construct, Row, panel, revcomp, strand_row

ROLES = ("Read 1", "Index 1 (i7)", "Index 2 (i5)", "Read 2")


@dataclass(frozen=True)
class SeqPrimer:
    role: str                    # one of ROLES
    name: str
    seq: str                     # a reference to a lib constant, never a literal in a page
    source: str = ""
    note: str = ""
    expect_mismatch: str = ""    # non-empty: the primer is NOT expected to match exactly


def custom(role: str, name: str, seq: str, source: str, note: str = "") -> SeqPrimer:
    """A protocol-specific primer. `seq` must still be a constant from the protocol module."""
    return SeqPrimer(role, name, seq, source, note)


def mismatching(p: SeqPrimer, why: str) -> SeqPrimer:
    return replace(p, expect_mismatch=why)


_ILL = 'Illumina "Illumina Adapter Sequences" #1000000002694'
_IDX = 'Illumina "Indexed Sequencing Overview Guide" #15057455'

TRUSEQ = {
    "R1": SeqPrimer("Read 1", "TruSeq Read 1", il.TRUSEQ_READ1, _ILL),
    "I1": SeqPrimer("Index 1 (i7)", "TruSeq Index 1", il.INDEX1_PRIMER, _ILL),
    "I2": SeqPrimer("Index 2 (i5)", "TruSeq Index 2, reverse-complement workflow",
                    il.INDEX2_PRIMER_RC, _IDX,
                    "Forward-strand instruments prime i5 off the flow-cell P5 oligo instead."),
    "R2": SeqPrimer("Read 2", "TruSeq Read 2", il.TRUSEQ_READ2, _ILL),
}
NEXTERA = {
    "R1": SeqPrimer("Read 1", "Nextera Read 1", nx.READ1_PRIMER, _ILL),
    "I1": SeqPrimer("Index 1 (i7)", "Nextera Index 1", nx.INDEX1_PRIMER, _ILL),
    "I2": SeqPrimer("Index 2 (i5)", "Nextera Index 2, reverse-complement workflow",
                    nx.INDEX2_PRIMER, _IDX,
                    "Forward-strand instruments prime i5 off the flow-cell P5 oligo instead."),
    "R2": SeqPrimer("Read 2", "Nextera Read 2", nx.READ2_PRIMER, _ILL),
}


@dataclass(frozen=True)
class Landing:
    strand: str        # "bottom" (primer = top-strand sequence, anneals to bottom) or "top"
    start: int         # primer site on the top-strand coordinate, half-open
    end: int
    covers: tuple[str, ...]
    reads: str         # first bases reported, 5'->3' of the read
    reads_from: str    # segment the read starts in
    free5: int = 0     # 5'-terminal primer bases with no exact partner (harmless: extension
                       # starts at the 3' end, which must match)


def _seg_at(lib: Construct, pos: int) -> str:
    p = 0
    for s in lib:
        if p <= pos < p + len(s):
            return s.name or "-"
        p += len(s)
    return "(off the end)"


def _covers(lib: Construct, a: int, b: int) -> tuple[str, ...]:
    out, p = [], 0
    for s in lib:
        if p < b and p + len(s) > a and s.name and s.name not in out:
            out.append(s.name)
        p += len(s)
    return tuple(out)


MIN_ANNEAL = 18    # a primer must pair over at least this many 3'-terminal bases; any
                   # 5'-terminal remainder is a flap (reported, harmless for extension)


def _hits(top: str, seq: str) -> list[tuple[str, int, int]]:
    hits = []
    i = top.find(seq)
    while i >= 0:
        hits.append(("bottom", i, i + len(seq)))
        i = top.find(seq, i + 1)
    rc = revcomp(seq)
    i = top.find(rc)
    while i >= 0:
        hits.append(("top", i, i + len(rc)))
        i = top.find(rc, i + 1)
    return hits


def locate_all(lib: Construct, p: SeqPrimer, n: int = 12) -> list[Landing]:
    """Every best-length site of ``p`` on ``lib``.

    Repeated sites are intentional in rolling-circle products. Ordinary sequencing
    libraries should call :func:`locate`, which continues to require uniqueness.
    """
    top = lib.top()
    hits = []
    for free5 in range(0, max(1, len(p.seq) - MIN_ANNEAL + 1)):
        hits = _hits(top, p.seq[free5:])
        if hits:
            break
    out = []
    for strand, a, b in hits:
        if strand == "bottom":   # same sense as top: extends rightwards, reports top
            reads, frm = top[b:b + n], _seg_at(lib, b)
        else:                    # anneals to top: extends leftwards, reports bottom 5'->3'
            # from the drawn bottom strand, so placeholders stay stand-ins, not bases
            reads, frm = lib.bottom()[max(0, a - n):a][::-1], _seg_at(lib, a - 1)
        out.append(Landing(strand, a, b, _covers(lib, a, b), reads, frm, free5))
    return out


def locate(lib: Construct, p: SeqPrimer, n: int = 12) -> Landing | None:
    """Unique site of ``p`` on the library, or ``None``."""
    hits = locate_all(lib, p, n)
    if len(hits) > 1:
        raise ValueError(f"{p.name}: {len(hits)} exact sites on {lib.name!r} -- ambiguous")
    return hits[0] if hits else None


def verify(lib: Construct, primers, required_roles=ROLES) -> list[str]:
    """Problems with declared primer sites and required run roles (empty = all good).

    `required_roles` describes the actual run. Single-end or single-index methods pass
    only the roles they use instead of inventing primers for absent reads.
    """
    errs = []
    for p in primers:
        hit = locate(lib, p)
        if hit is None and not p.expect_mismatch:
            errs.append(f"{p.role} / {p.name}: no exact site on {lib.name!r}")
        if hit is not None and p.expect_mismatch:
            errs.append(f"{p.role} / {p.name}: declared mismatching but matches exactly")
    roles = {p.role for p in primers}
    missing = [r for r in required_roles if r not in roles]
    if missing:
        errs.append(f"no primer declared for {', '.join(missing)}")
    return errs


def diagram(lib: Construct, primers, caption: str = "Sequencing primers on the final library") -> str:
    """Draw every declared sequencing primer at its computed site on the final duplex."""
    located = [(p, locate(lib, p)) for p in primers]
    gutter = len("final library  ")

    def primer_row(p: SeqPrimer, hit: Landing) -> Row:
        role = f"{p.role} / {p.name}"
        tag = "r1" if p.role == "Read 1" else "r2" if p.role == "Read 2" else "cbc"
        if hit.strand == "bottom":
            # Same sense as the top strand; a 5' flap extends to the left of the match.
            start = hit.start - hit.free5
            return Row(chunks=[(p.seq, tag, False)], indent=gutter + start, prefix="5'- ",
                       suffix=f" -3'  --------> {role}")
        # Antiparallel below the top strand.  Reverse for the left-to-right 3'->5' drawing;
        # a 5' flap consequently extends to the right of the exact match.
        return Row(chunks=[(p.seq[::-1], tag, False)], indent=gutter + hit.start, prefix="3'- ",
                   suffix=f" -5'  <-------- {role}")

    above = [primer_row(p, h) for p, h in located if h and h.strand == "bottom"]
    below = [primer_row(p, h) for p, h in located if h and h.strand == "top"]
    missing = [Row(chunks=[(f"{p.role} / {p.name}: no exact site — {p.expect_mismatch}",
                            None, False)]) for p, h in located if h is None]
    top = strand_row(lib, "top", prefix="final library  5'- ", suffix=" -3'")
    bottom = strand_row(lib, "bottom", prefix="               3'- ", suffix=" -5'")
    return panel([*above, top, bottom, *below, *missing], cls="long", caption=caption)


def unavailable_diagram(lib: Construct, reason: str,
                        roles: tuple[str, ...] = ("Read 1", "Index", "Read 2"),
                        caption: str = "Sequencing-primer binding on the final library") -> str:
    """Draw the final duplex and make an unsupported primer placement visibly absent."""
    rows = [strand_row(lib, "top", prefix="final library  5'- ", suffix=" -3'"),
            strand_row(lib, "bottom", prefix="               3'- ", suffix=" -5'")]
    rows.extend(Row(chunks=[(f"{role} primer: ?  {reason}", None, False)]) for role in roles)
    return panel(rows, cls="long", caption=caption)


def section(lib: Construct, primers, heading: str = "Sequencing primers",
            intro: str = "", required_roles=ROLES) -> str:
    """Render the page section. Raises if `verify` finds anything."""
    errs = verify(lib, primers, required_roles=required_roles)
    if errs:
        raise ValueError("sequencing primers: " + "; ".join(errs))
    order = {r: i for i, r in enumerate(ROLES)}
    rows = []
    for p in sorted(primers, key=lambda p: order[p.role]):
        hit = locate(lib, p)
        if hit:
            where = (f"anneals to the {hit.strand} strand over "
                     + ", ".join(html.escape(c) for c in hit.covers))
            if hit.free5:
                where += (f"; its 5'-terminal {hit.free5} nt are an unpaired flap "
                          "(harmless: extension starts at the paired 3' end)")
            first = (f"<code>{html.escape(hit.reads)}</code>&hellip; "
                     f"(from {html.escape(hit.reads_from)})")
        else:
            where = f"<b>no exact site</b> &mdash; {html.escape(p.expect_mismatch)}"
            first = "&mdash;"
        note = f"<br><small>{html.escape(p.note)}</small>" if p.note else ""
        src = f"<br><small>{html.escape(p.source)}</small>" if p.source else ""
        rows.append(f"<tr><td>{html.escape(p.role)}</td><td>{html.escape(p.name)}{src}</td>"
                    f"<td><code>5'-{p.seq}-3'</code> ({len(p.seq)} nt){note}</td>"
                    f"<td>{where}</td><td>{first}</td></tr>")
    intro_html = f"<p><info>{intro}</info></p>\n" if intro else ""
    return (f"<h2>{html.escape(heading)}</h2>\n{intro_html}"
            "<p><info>Sequences come from <code>lib/</code>; the landing site and the first "
            "bases of each read are computed from this page's final library, so they cannot "
            "drift from it.</info></p>\n" + diagram(lib, primers) + "\n<div class=\"tw\"><table>\n"
            "<tr><th>Read</th><th>Primer</th><th>Sequence</th><th>Lands on</th>"
            "<th>First bases read</th></tr>\n" + "\n".join(rows) + "\n</table></div>\n")
