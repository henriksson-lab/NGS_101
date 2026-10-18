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
from chemdraw import (ArrowVisual, Construct, MolecularFeature, Row, Segment, StrandVisual,
                      feature_rows, inline_sequence, panel, revcomp, strand_row)

ROLES = ("Read 1", "Index 1 (i7)", "Index 2 (i5)", "Read 2")
READ_USES = {
    "include": "Yes — include in sequence analysis",
    "extract": "No — extract identifiers only",
    "extract_then_include": "Yes — after extracting identifiers",
    "demultiplex": "No — use for sample demultiplexing",
    "extract_and_demultiplex": "No — extract identifiers and demultiplex",
    "technical": "No — technical sequence only",
    "not_used": "No — acquired but not used for sequence analysis",
    "not_applicable": "No — assay readout is not sequencing",
    "undetermined": "Undetermined",
}


@dataclass(frozen=True)
class SeqPrimer:
    role: str                    # one of ROLES
    name: str
    seq: str                     # a reference to a lib constant, never a literal in a page
    source: str = ""
    note: str = ""
    expect_mismatch: str = ""    # non-empty: the primer is NOT expected to match exactly
    min_anneal: int | None = None # override only for a source-established shorter 3′ site


def custom(role: str, name: str, seq: str, source: str, note: str = "",
           *, min_anneal: int | None = None) -> SeqPrimer:
    """A protocol-specific primer. `seq` must still be a constant from the protocol module."""
    return SeqPrimer(role, name, seq, source, note, min_anneal=min_anneal)


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


def standard_primers(lib: Construct, acquired_roles=ROLES) -> tuple[SeqPrimer, ...]:
    """Infer the standard Illumina primers used by an acquired run from ``lib``.

    A protocol supplies only which reads were acquired.  Primer family, binding strand
    and direction come from sequence matching against the completed construct.  If no
    standard primer—or more than one distinct standard primer—can supply a requested
    role, construction fails and the protocol must declare its custom primer explicitly.
    """
    roles = tuple(acquired_roles)
    unknown = set(roles) - set(ROLES)
    if unknown:
        raise ValueError(f"unknown acquired read roles: {', '.join(sorted(unknown))}")
    out = []
    standards = (*TRUSEQ.values(), *NEXTERA.values())
    for role in roles:
        by_sequence = {}
        for primer in standards:
            if primer.role == role and locate(lib, primer) is not None:
                by_sequence.setdefault(primer.seq, primer)
        if not by_sequence:
            raise ValueError(f"no standard Illumina primer binds for {role}")
        if len(by_sequence) != 1:
            names = ", ".join(p.name for p in by_sequence.values())
            raise ValueError(f"multiple standard Illumina primers bind for {role}: {names}")
        out.append(next(iter(by_sequence.values())))
    return tuple(out)


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


@dataclass(frozen=True)
class FeatureReadSpan:
    """One annotated molecular feature observed in a declared sequencing read."""
    read: str
    feature: MolecularFeature
    segment: str
    cycle_start: int
    cycle_end: int

    @property
    def cycles(self) -> str:
        return (str(self.cycle_start) if self.cycle_start == self.cycle_end else
                f"{self.cycle_start}–{self.cycle_end}")


@dataclass(frozen=True)
class ReadContent:
    """One run read, with content and handling derived from library geometry."""
    read: str
    cycles: str
    content: str
    handling: str


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


def _read_preview(lib: Construct, start: int, strand: str, n: int) -> str:
    """Show known first bases, stopping at a schematic region of unknown bases."""
    end = min(len(lib), start + n) if strand == "bottom" else max(0, start - n)
    position = 0
    unknown = False
    for segment in lib:
        after = position + len(segment)
        if segment.is_role_token():
            if strand == "bottom" and start <= position < end:
                end, unknown = position, True
                break
            if strand == "top" and end < after <= start:
                end, unknown = after, True
                break
        position = after
    if strand == "bottom":
        bases = lib.top()[start:end]
    else:
        bases = lib.bottom()[end:start][::-1]
    return bases + ("…" if unknown else "")


def locate_all(lib: Construct, p: SeqPrimer, n: int = 12,
               min_anneal: int = MIN_ANNEAL) -> list[Landing]:
    """Every best-length site of ``p`` on ``lib``.

    Repeated sites are intentional in rolling-circle products. Ordinary sequencing
    libraries should call :func:`locate`, which continues to require uniqueness.
    """
    top = lib.top()
    hits = []
    if p.min_anneal is not None:
        min_anneal = p.min_anneal
    if min_anneal < 1:
        raise ValueError("minimum annealing length must be positive")
    # A primer shorter than the general threshold can still have one valid site: its
    # entire sequence must match.  This preserves exact placement of short custom
    # sequencing primers while preventing any 5′ trimming below the requested footprint.
    effective_min = min(min_anneal, len(p.seq))
    for free5 in range(0, len(p.seq) - effective_min + 1):
        hits = _hits(top, p.seq[free5:])
        if hits:
            break
    out = []
    for strand, a, b in hits:
        if strand == "bottom":   # same sense as top: extends rightwards, reports top
            reads, frm = _read_preview(lib, b, strand, n), _seg_at(lib, b)
        else:                    # anneals to top: extends leftwards, reports bottom 5'->3'
            reads, frm = _read_preview(lib, a, strand, n), _seg_at(lib, a - 1)
        out.append(Landing(strand, a, b, _covers(lib, a, b), reads, frm, free5))
    return out


def locate(lib: Construct, p: SeqPrimer, n: int = 12,
           min_anneal: int = MIN_ANNEAL) -> Landing | None:
    """Unique site of ``p`` on the library, or ``None``."""
    hits = locate_all(lib, p, n, min_anneal)
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


def _cycle_coordinates(lib: Construct, hit: Landing):
    """Convert drawing coordinates to molecular coordinates for read-cycle arithmetic.

    Known placeholder lengths override drawing widths. Unknown lengths retain temporary
    coordinates only so callers can stop at them; callers must never report spans across
    them. A primer boundary inside a compressed placeholder is not an exact molecular site.
    """
    coords, boundaries = [], {}
    drawing = molecular = 0
    for segment in lib:
        width = len(segment)
        length = segment.length_bp if segment.length_bp is not None else width
        for boundary in (hit.start, hit.end):
            if drawing <= boundary <= drawing + width:
                offset = boundary - drawing
                if offset not in (0, width) and length != width:
                    raise ValueError("sequencing primer boundary lies inside a compressed placeholder")
                boundaries[boundary] = molecular + (length if offset == width else offset)
        coords.append((segment, molecular, molecular + length))
        drawing += width
        molecular += length
    return coords, replace(hit, start=boundaries[hit.start], end=boundaries[hit.end])


def feature_spans(lib: Construct, primers, read_lengths: dict[str, int]) -> list[FeatureReadSpan]:
    """Map annotated regions to read cycles from primer geometry and declared run lengths.

    Cycle coordinates are never stored on a segment: they follow from the final construct,
    the primer landing strand, and the protocol's actual number of cycles.  Features beyond
    a read's declared length are omitted; a feature crossing the end is clipped.
    """
    unknown = set(read_lengths) - {p.role for p in primers}
    if unknown:
        raise ValueError(f"read lengths name undeclared roles: {', '.join(sorted(unknown))}")
    if any(not isinstance(n, int) or n <= 0 for n in read_lengths.values()):
        raise ValueError("read lengths must be positive integers")
    def crosses_unknown(coords, a: int, b: int) -> bool:
        """Whether a cycle path crosses a schematic region of unknown real length.

        Role tokens have a printable width so diagrams can show them, but that width is
        not a molecular length.  Identifier cycles beyond such a token are unknowable and
        must not be derived from the number of characters in its label.
        """
        lo, hi = sorted((a, b))
        return any(_unstated_length(s) and start < hi and end > lo
                   for s, start, end in coords)
    out = []
    for primer in primers:
        limit = read_lengths.get(primer.role)
        if limit is None:
            continue
        hit = locate(lib, primer)
        if hit is None:
            continue
        coords, hit = _cycle_coordinates(lib, hit)
        for segment, start, end in coords:
            if not segment.feature:
                continue
            if hit.strand == "bottom":
                if end <= hit.end:
                    continue
                if crosses_unknown(coords, hit.end, start) or (
                        segment.is_role_token() and segment.length_bp is None):
                    continue
                first = max(start, hit.end) - hit.end + 1
                last = end - hit.end
            else:
                if start >= hit.start:
                    continue
                if crosses_unknown(coords, end, hit.start) or (
                        segment.is_role_token() and segment.length_bp is None):
                    continue
                first = hit.start - min(end, hit.start) + 1
                last = hit.start - start
            if first > limit:
                continue
            out.append(FeatureReadSpan(primer.role, segment.feature, segment.name,
                                       first, min(last, limit)))
    return sorted(out, key=lambda x: (ROLES.index(x.read), x.cycle_start, x.feature.id))


def _unstated_length(segment: Segment) -> bool:
    """Whether a placeholder's drawing width is not a molecular length."""
    return (segment.placeholder and segment.length_bp is None
            and (not segment.feature or segment.is_role_token()))


UNPLACED_SUFFIX = "acquired; priming site not represented on this library"
UNSPECIFIED_HANDLING = "Undetermined — acquisition not stated"


def _unplaced_by_role(unplaced_reads, primer_roles) -> dict[str, ReadContent]:
    """Validate source-declared acquired reads whose priming site is not represented."""
    out = {}
    for row in unplaced_reads or ():
        if not isinstance(row, ReadContent):
            raise TypeError("an unplaced read must be a read_spec() row")
        if row.read not in ROLES:
            raise ValueError(f"unplaced read needs a run read role, not {row.read!r}")
        if row.read in out:
            raise ValueError(f"unplaced read {row.read!r} is declared twice")
        if row.read in primer_roles:
            raise ValueError(f"{row.read} has a declared primer; it is not unplaced")
        out[row.read] = row
    return out


def read_contents(lib: Construct, primers,
                  read_lengths: dict[str, int] | None = None,
                  read_uses: dict[str, str] | None = None,
                  acquired_roles=None, *, unplaced_reads=(),
                  acquisition_unspecified: bool = False) -> list[ReadContent]:
    """Describe R1/I1/I2/R2 from the construct and declared primer sites.

    The construct supplies order and meaning; primer geometry supplies direction.  A
    protocol need only declare actual run lengths when the source fixes them.  Missing
    primer roles and declared primer sites omitted from ``acquired_roles`` become explicit
    ``not acquired`` rows rather than disappearing.

    ``unplaced_reads`` are :func:`read_spec` rows for reads the source says were acquired
    although their priming site is not represented (for example an unprinted custom
    primer).  They render as acquired without a site, never as ``not acquired``.
    ``acquisition_unspecified`` keeps every declared primer site and its derived content
    but states that the source does not say which reads were acquired.
    """
    primers = tuple(primers)
    read_lengths = read_lengths or {}
    read_uses = read_uses or {}
    unknown = set(read_lengths) - {p.role for p in primers}
    if unknown:
        raise ValueError(f"read lengths name undeclared roles: {', '.join(sorted(unknown))}")
    unknown_uses = set(read_uses) - {p.role for p in primers}
    if unknown_uses:
        raise ValueError(f"read uses name undeclared roles: {', '.join(sorted(unknown_uses))}")
    invalid_uses = set(read_uses.values()) - set(READ_USES)
    if invalid_uses:
        raise ValueError(f"unknown read use: {', '.join(sorted(invalid_uses))}")
    by_role = {p.role: p for p in primers}
    if len(by_role) != len(primers):
        raise ValueError("at most one sequencing primer may be declared for each read role")
    unplaced = _unplaced_by_role(unplaced_reads, set(by_role))
    if acquisition_unspecified and (acquired_roles or unplaced):
        raise ValueError("unspecified acquisition cannot also declare acquired reads")
    acquired = set(by_role) if acquired_roles is None else set(acquired_roles)
    acquired -= set(unplaced)
    if acquisition_unspecified:
        acquired = set(by_role)
    undeclared_acquired = acquired - set(by_role)
    if undeclared_acquired:
        raise ValueError("acquired reads lack a declared primer: "
                         + ", ".join(sorted(undeclared_acquired)))

    out = []
    for role in ROLES:
        primer = by_role.get(role)
        if role in unplaced:
            row = unplaced[role]
            out.append(replace(row, content=f"{row.content} ({UNPLACED_SUFFIX})"))
            continue
        if primer is None and acquisition_unspecified:
            out.append(ReadContent(role, "—", "no sequencing primer declared",
                                   UNSPECIFIED_HANDLING))
            continue
        if primer is None or role not in acquired:
            out.append(ReadContent(role, "—", "not acquired", "No — not sequenced"))
            continue
        hit = locate(lib, primer)
        if hit is None:
            out.append(ReadContent(role, "not fixed", "cannot be derived from an exact site",
                                   "Undetermined"))
            continue
        coords, hit = _cycle_coordinates(lib, hit)
        limit = read_lengths.get(role)
        candidates = []
        for segment, start, end in coords:
            if hit.strand == "bottom":
                if end <= hit.end:
                    continue
                first = max(start, hit.end) - hit.end + 1
                last = end - hit.end
            else:
                if start >= hit.start:
                    continue
                first = hit.start - min(end, hit.start) + 1
                last = hit.start - start
            if limit is not None and first > limit:
                continue
            last = min(last, limit) if limit is not None else last
            if first > last:
                continue
            candidates.append((segment, first, last))

        # Reading order is reverse genomic order for a primer annealed to the top strand.
        candidates.sort(key=lambda x: x[1])
        path = []
        following_unresolved = None
        short_molecule_readthrough = False
        for candidate_no, (segment, first, last) in enumerate(candidates):
            # Identifier placeholders such as B*16 have a real, modelled length even
            # though their bases vary.  An unfeatured payload placeholder such as
            # ``cDNA insert`` is the molecule-dependent case.  A featured role token
            # (for example ``[i7 index; length not stated]``) has drawing width but no
            # stated length, so later cycles are unknowable, as in feature_spans().
            unknown_length = _unstated_length(segment)
            path.append((segment, first, last))
            if unknown_length:
                more = candidate_no + 1 < len(candidates)
                short_molecule_readthrough = more and not segment.feature
                index_stop = (limit is None and role.startswith("Index") and segment.feature
                              and segment.feature.role == "sample_index")
                if more and (segment.inferred or segment.feature) and not index_stop:
                    following_unresolved = candidates[candidate_no + 1][0]
                break
            if (limit is None and role.startswith("Index")
                                   and segment.feature
                                   and segment.feature.role == "sample_index"):
                break
        parts = []
        for segment, first, last in path:
            label = segment.feature.label if segment.feature else (segment.name or "unnamed region")
            span = str(first) if first == last else f"{first}–{last}"
            parts.append(f"{label} (cycles {span})")
            # A schematic payload or role-token width must not expose later adapter
            # sequence or imply a cycle count merely because the placeholder is drawn.
            if _unstated_length(segment):
                kind = "length not stated" if segment.feature else "molecule-dependent length"
                parts[-1] = f"{label} (from cycle {first}; {kind})"
                break
        if short_molecule_readthrough:
            parts.append("short molecules may read through into the opposite library end")
        if following_unresolved is not None:
            label = (following_unresolved.feature.label if following_unresolved.feature
                     else (following_unresolved.name or "unnamed region"))
            parts.append(f"then {label} (cycle boundary unresolved)")

        features = {segment.feature.role for segment, _, _ in path if segment.feature}
        identifier = bool(features - {"sample_index"})
        technical_tags = {"p5", "p7", "s5", "s7", "me", "r1", "r2", "tso",
                          "cbc", "umi", "insert"}
        payload = any(not segment.feature and not set((segment.tag or "").split("+"))
                      & technical_tags for segment, _, _ in path)
        if role.startswith("Index") and identifier and "sample_index" in features:
            handling = "No — extract identifiers and demultiplex"
        elif role.startswith("Index") and identifier:
            handling = "No — extract identifiers only"
        elif role.startswith("Index") and "sample_index" in features:
            handling = "No — use for sample demultiplexing"
        elif role.startswith("Index") and payload:
            handling = "Undetermined — non-identifier content in an index read"
        elif role.startswith("Index"):
            handling = "No — technical sequence only"
        elif identifier and payload:
            handling = "Yes — after extracting identifiers"
        elif identifier:
            handling = "No — extract identifiers only"
        elif payload:
            handling = "Yes — include in sequence analysis"
        else:
            handling = "No — technical sequence only"
        if role in read_uses:
            handling = READ_USES[read_uses[role]]
        if acquisition_unspecified:
            handling = UNSPECIFIED_HANDLING
        cycles = str(limit) if limit is not None else "not fixed"
        out.append(ReadContent(role, cycles, "; ".join(parts) or "read-through beyond library",
                               handling))
    return out


def read_content_table(lib: Construct, primers,
                       read_lengths: dict[str, int] | None = None,
                       read_uses: dict[str, str] | None = None,
                       acquired_roles=None, *, unplaced_reads=(),
                       acquisition_unspecified: bool = False) -> str:
    """Render the four standard run reads, including roles that are not acquired."""
    rows = read_contents(lib, tuple(primers), read_lengths, read_uses, acquired_roles,
                         unplaced_reads=unplaced_reads,
                         acquisition_unspecified=acquisition_unspecified)
    body = "\n".join(
        f"<tr><td>{html.escape(x.read)}</td><td>{html.escape(x.cycles)}</td>"
        f"<td>{html.escape(x.content)}</td><td>{html.escape(x.handling)}</td></tr>"
        for x in rows)
    extra = ""
    if unplaced_reads:
        extra += (" Rows marked as acquired without a represented priming site are "
                  "source-declared, not derived from the drawing.")
    if acquisition_unspecified:
        extra += (" The source does not state which reads were acquired; each row shows "
                  "what a declared primer site would read.")
    return ("<h3>Read content and handling</h3>\n"
            "<p><info>Content and cycle positions are derived from the final construct, "
            "primer direction, molecular-feature annotations and declared run lengths. "
            "\u2018No\u2019 means the read is auxiliary or absent, not that its information is "
            "unimportant." + extra + "</info></p>\n<div class=\"tw\"><table>\n"
            "<tr><th>Read</th><th>Cycles</th><th>Content</th>"
            "<th>Include in sequence analysis?</th></tr>\n" + body + "\n</table></div>\n")


def read_spec(read: str, cycles: int | str | None, content: str,
              use: str) -> ReadContent:
    """Construct one source-declared read row using the controlled handling vocabulary."""
    if read not in ROLES and read != "Sequencing output":
        raise ValueError(f"unknown read role: {read}")
    if use not in READ_USES:
        raise ValueError(f"unknown read use: {use}")
    if isinstance(cycles, int) and cycles <= 0:
        raise ValueError("read cycles must be positive")
    cycle_text = "not fixed" if cycles is None else str(cycles)
    return ReadContent(read, cycle_text, content, READ_USES[use])


PARTIAL_GEOMETRY_NOTE = ("Read roles and contents are source-declared; computed primer "
                         "landing is shown above only for the primer sites this model "
                         "represents.")


def read_spec_table(rows, note: str = "") -> str:
    """Render structured source-declared reads when exact primer geometry is unavailable.

    The default note claims that no exact primer geometry is available.  A caller that
    has drawn some primer sites passes ``PARTIAL_GEOMETRY_NOTE`` or its own wording.
    """
    rows = tuple(rows)
    if not rows:
        raise ValueError("a read-content table needs at least one row")
    body = "\n".join(
        f"<tr><td>{html.escape(x.read)}</td><td>{html.escape(x.cycles)}</td>"
        f"<td>{html.escape(x.content)}</td><td>{html.escape(x.handling)}</td></tr>"
        for x in rows)
    explanation = (note or "Read roles and contents are source-declared because exact "
                   "sequencing-primer geometry is unavailable.")
    return ("<h3>Read content and handling</h3>\n"
            f"<p><info>{html.escape(explanation)}</info></p>\n<div class=\"tw\"><table>\n"
            "<tr><th>Read</th><th>Cycles</th><th>Content</th>"
            "<th>Include in sequence analysis?</th></tr>\n" + body + "\n</table></div>\n")


def nanopore_read_content_table(lib: Construct | None,
                                adapter_name: str = "motor adapter") -> str:
    """Describe primer-free nanopore output without inventing a primer site.

    Nanopore libraries are metered through a pore by a motor-bearing adapter; treating
    the absent synthesis primer as an unresolved oligo is a category error.  The exact
    motor/adaptor bases may still be proprietary, but that uncertainty belongs to the
    final construct rather than to read initiation.
    """
    labels = []
    payloads = []
    technical_tags = {"p5", "p7", "s5", "s7", "me", "r1", "r2", "r3", "tso",
                      "cbc", "umi", "insert", "t7", "w1"}
    for segment in lib or ():
        label = segment.feature.label if segment.feature else ""
        if label and label not in labels:
            labels.append(label)
        tags = set((segment.tag or "").split("+"))
        if not segment.feature and not tags & technical_tags and segment.name:
            if segment.name not in payloads and "adapter" not in segment.name.lower():
                payloads.append(segment.name)
    content = ", ".join([*labels, *payloads[:4]]) or "the library strand"
    if len(payloads) > 4:
        content += ", other library regions"
    return read_spec_table(
        [read_spec("Sequencing output", None,
                   f"One strand enters from the {adapter_name} end and reports {content}",
                   "include")],
        "Nanopore sequencing is motor-driven and does not use a synthesis sequencing primer.")


def ensure_read_content(page_html: str, reason: str) -> str:
    """Ensure a generated protocol page cannot silently omit its read description.

    Protocol builders should provide geometry-derived or source-declared rows.  This is
    the final publication boundary: it adds an explicit unresolved row when a legacy or
    source-dependent builder provides neither.  It deliberately does not infer read names,
    orientation or handling from prose.
    """
    if "<h3>Read content and handling</h3>" in page_html:
        return page_html
    fallback = read_spec_table(
        [read_spec("Sequencing output", None,
                   "Read structure is not yet represented as structured data on this page",
                   "undetermined")],
        "Structured read placement is unresolved: " + reason,
    )
    closing = page_html.rfind("</div>")
    if closing < 0:
        return page_html + "\n" + fallback
    return page_html[:closing] + fallback + "\n" + page_html[closing:]


def unresolved_read_content_table(lib: Construct | None, reason: str,
                                  read_lengths: dict[str, int] | None = None,
                                  acquired_roles=None) -> str:
    """Honest fallback when exact read-entry geometry is unavailable.

    When a final construct exists, content is summarized from structured molecular
    features and non-technical regions.  Platform workflows whose proprietary adapter
    sequence prevents a final construct may pass ``None``; they still receive an explicit
    generic sequencing-output row.  No cycle placement or read direction is invented.
    Protocols should replace this with :func:`read_spec_table` once a source establishes
    their run layout.
    """
    technical_tags = {"p5", "p7", "s5", "s7", "me", "r1", "r2", "r3", "tso",
                      "cbc", "umi", "insert", "t7", "w1"}
    labels = []
    for segment in lib or ():
        label = segment.feature.label if segment.feature else ""
        if label and label not in labels:
            labels.append(label)
    payloads = []
    for segment in lib or ():
        tags = set((segment.tag or "").split("+"))
        if segment.feature or tags & technical_tags or not segment.name:
            continue
        if segment.name not in payloads:
            payloads.append(segment.name)
    shown_payloads = payloads[:4]
    if len(payloads) > len(shown_payloads):
        shown_payloads.append("other library regions")
    summary = (", ".join([*labels, *shown_payloads]) or
               ("final-library molecule" if lib is not None else "sequencing output"))
    summary += " (read orientation and cycle boundaries unresolved)"
    if labels and payloads:
        use = "extract_then_include"
    elif labels:
        use = "extract"
    elif payloads:
        use = "include"
    else:
        use = "technical"
    read_lengths = read_lengths or {}
    acquired = tuple(acquired_roles or read_lengths)
    if acquired:
        unknown = set(acquired) - set(ROLES)
        if unknown:
            raise ValueError(f"unknown acquired read roles: {', '.join(sorted(unknown))}")
        rows = [read_spec(role, read_lengths.get(role), summary, use) for role in ROLES
                if role in acquired]
        rows.extend(ReadContent(role, "—", "not acquired", "No — not sequenced")
                    for role in ROLES if role not in acquired)
        rows.sort(key=lambda x: ROLES.index(x.read))
    else:
        rows = [read_spec("Sequencing output", None, summary, use)]
    return read_spec_table(rows, "Exact read placement is unresolved: " + reason)


def diagram(lib: Construct, primers, caption: str = "Sequencing primers on the final library") -> str:
    """Draw every declared sequencing primer at its computed site on the final duplex."""
    located = [(p, locate(lib, p)) for p in primers]
    label = "final library  "
    sequence_column = len(label + "5'- ")

    def primer_row(p: SeqPrimer, hit: Landing) -> Row:
        role = f"{p.role} / {p.name}"
        tag = "r1" if p.role == "Read 1" else "r2" if p.role == "Read 2" else "cbc"
        if hit.strand == "bottom":
            # Same sense as the top strand; a 5' flap extends to the left of the match.
            start = hit.start - hit.free5
            primer = Construct([Segment(p.name, p.seq, tag)], name=p.name)
            row = strand_row(primer, "top", indent=len(label) + start,
                             prefix="5'- ", suffix=" -3'")
            strand = row.visual
            assert isinstance(strand, StrandVisual)
            arrow_start = strand.end + len(row.suffix) + 2
            row.visual = (strand, ArrowVisual(arrow_start, arrow_start + 8, "right", role))
            return row
        # Antiparallel below the top strand.  Reverse for the left-to-right 3'->5' drawing;
        # a 5' flap consequently extends to the right of the exact match.
        primer = Construct([Segment(p.name, revcomp(p.seq), tag)], name=p.name)
        row = strand_row(primer, "bottom", indent=len(label) + hit.start,
                         prefix="3'- ", suffix=" -5'")
        strand = row.visual
        assert isinstance(strand, StrandVisual)
        arrow_start = strand.end + len(row.suffix) + 2
        row.visual = (strand, ArrowVisual(arrow_start, arrow_start + 8, "left", role,
                                          label_after=True))
        return row

    above = [primer_row(p, h) for p, h in located if h and h.strand == "bottom"]
    below = [primer_row(p, h) for p, h in located if h and h.strand == "top"]
    missing = [Row(chunks=[(f"{p.role} / {p.name}: no exact site — {p.expect_mismatch}",
                            None, False)]) for p, h in located if h is None]
    top = strand_row(lib, "top", prefix="final library  5'- ", suffix=" -3'")
    bottom = strand_row(lib, "bottom", prefix="               3'- ", suffix=" -5'")
    semantics = feature_rows(lib, prefix_width=sequence_column)
    return panel([*above, top, bottom, *below, *semantics, *missing],
                 cls="long", caption=caption)


def unavailable_diagram(lib: Construct, reason: str,
                        roles: tuple[str, ...] = ("Read 1", "Index", "Read 2"),
                        caption: str = "Sequencing-primer binding on the final library") -> str:
    """Draw the final duplex and make an unsupported primer placement visibly absent."""
    rows = [strand_row(lib, "top", prefix="final library  5'- ", suffix=" -3'"),
            strand_row(lib, "bottom", prefix="               3'- ", suffix=" -5'")]
    rows.extend(feature_rows(lib, prefix_width=len("final library  5'- ")))
    rows.extend(Row(chunks=[(f"{role} primer: ?  {reason}", None, False)]) for role in roles)
    return panel(rows, cls="long", caption=caption)


def section(lib: Construct, primers, heading: str = "Sequencing primers",
            intro: str = "", required_roles=ROLES,
            read_lengths: dict[str, int] | None = None,
            read_uses: dict[str, str] | None = None,
            content_primers=None,
            show_read_content: bool = True, *, unplaced_reads=(),
            acquisition_unspecified: bool = False) -> str:
    """Render primers and reads; ``required_roles`` is also the set actually acquired.

    ``content_primers`` selects the loaded primer when ``primers`` deliberately includes
    alternatives or counterexamples for explanation.  The full set remains in the primer
    table, while read content is derived only from the selected run primers.
    """
    unplaced_roles = {row.read for row in unplaced_reads or ()}
    required_roles = (() if acquisition_unspecified else
                      tuple(r for r in required_roles if r not in unplaced_roles))
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
            continuation = "" if hit.reads.endswith("…") else "&hellip;"
            first = (f"<code>{html.escape(hit.reads)}</code>{continuation} "
                     f"(from {html.escape(hit.reads_from)})")
        else:
            where = f"<b>no exact site</b> &mdash; {html.escape(p.expect_mismatch)}"
            first = "&mdash;"
        note = f"<br><small>{html.escape(p.note)}</small>" if p.note else ""
        src = f"<br><small>{html.escape(p.source)}</small>" if p.source else ""
        sequence = inline_sequence(p.seq, p.name)
        rows.append(f"<tr><td>{html.escape(p.role)}</td><td>{html.escape(p.name)}{src}</td>"
                    f"<td>{sequence}{note}</td>"
                    f"<td>{where}</td><td>{first}</td></tr>")
    intro_html = f"<p><info>{intro}</info></p>\n" if intro else ""
    feature_table = ""
    if read_lengths and not show_read_content:
        spans = feature_spans(lib, primers, read_lengths)
        if spans:
            body = "\n".join(
                f"<tr><td>{html.escape(x.read)}</td><td>{x.cycles}</td>"
                f"<td>{html.escape(x.feature.label)}</td>"
                f"<td>{html.escape(x.segment)}</td>"
                f"<td>{html.escape(x.feature.encoding)}</td></tr>" for x in spans
            )
            feature_table = ("\n<h3>Identifier cycles</h3>\n"
                "<p><info>Cycle ranges are computed from the final construct, primer "
                "landing direction and declared run length.</info></p>\n"
                "<div class=\"tw\"><table><tr><th>Read</th><th>Cycles</th>"
                "<th>Role</th><th>Region</th><th>Encoding</th></tr>\n" + body
                + "\n</table></div>\n")
    content_table = (read_content_table(lib, content_primers or primers,
                                        read_lengths, read_uses,
                                        acquired_roles=(None if acquisition_unspecified
                                                        else required_roles),
                                        unplaced_reads=unplaced_reads,
                                        acquisition_unspecified=acquisition_unspecified)
                     if show_read_content else "")
    return (f"<h2>{html.escape(heading)}</h2>\n{intro_html}"
            "<p><info>Sequences come from <code>lib/</code>; the landing site and the first "
            "bases of each read are computed from this page's final library, so they cannot "
            "drift from it.</info></p>\n" + diagram(lib, primers) + "\n<div class=\"tw\"><table>\n"
            "<tr><th>Read</th><th>Primer</th><th>Sequence</th><th>Lands on</th>"
            "<th>First bases read</th></tr>\n" + "\n".join(rows)
            + "\n</table></div>\n" + content_table + feature_table)
