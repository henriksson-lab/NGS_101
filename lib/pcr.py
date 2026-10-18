"""Checked primer-binding geometry for PCR reaction states.

A PCR product alone hides which ends were selected and where primer extension starts.
This module makes the binding state mandatory data for ``Workflow.pcr``. Exact sites are
located from sequence; region-level sites are reserved for published primer pools whose
target class is stated but whose transcript coordinates are not supplied by the source.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import seqprimers
from chemdraw import (ArrowVisual, Construct, Row, Scene, Segment, SpanVisual, StrandVisual,
                      revcomp, strand_row)


@dataclass(frozen=True)
class PCRPrimer:
    name: str
    sequence: str
    anneal_nt: int = seqprimers.MIN_ANNEAL

    def __post_init__(self) -> None:
        sequence = self.sequence.upper()
        if not sequence or set(sequence) - set("ACGTRYSWKMBDHVN"):
            raise ValueError(f"PCR primer {self.name!r} must be an IUPAC DNA sequence")
        if not 1 <= self.anneal_nt <= len(sequence):
            raise ValueError("PCR primer annealing length must fit within the primer")


def _landing(template: Construct, primer: PCRPrimer) -> seqprimers.Landing:
    query = seqprimers.SeqPrimer("PCR", primer.name, primer.sequence)
    hit = seqprimers.locate(template, query, min_anneal=primer.anneal_nt)
    if hit is None:
        raise ValueError(f"PCR primer {primer.name!r} has no site on {template.name!r}")
    return hit


def _landing_on(template: Construct, primer: PCRPrimer,
                strand: str) -> seqprimers.Landing:
    """Locate one primer species on one requested strand of a symmetric template."""
    query = seqprimers.SeqPrimer("PCR", primer.name, primer.sequence)
    hits = [hit for hit in seqprimers.locate_all(
        template, query, min_anneal=primer.anneal_nt) if hit.strand == strand]
    if len(hits) != 1:
        raise ValueError(
            f"PCR primer {primer.name!r} has {len(hits)} {strand}-strand sites "
            f"on {template.name!r}; expected exactly one")
    return hits[0]


def _primer_row(template: Construct, primer: PCRPrimer, hit: seqprimers.Landing,
                *, gutter: int) -> Row:
    tag = "r1" if hit.strand == "bottom" else "r2"
    if hit.strand == "bottom":
        shown = Construct([Segment(primer.name, primer.sequence, tag)], primer.name)
        row = strand_row(shown, "top", indent=gutter + hit.start - hit.free5,
                         prefix="5'- ", suffix=" -3'")
        strand = row.visual
        assert isinstance(strand, StrandVisual)
        row.visual = (strand, ArrowVisual(strand.end + len(row.suffix) + 2,
                                          strand.end + len(row.suffix) + 10,
                                          "right", primer.name))
        return row
    shown = Construct([Segment(primer.name, revcomp(primer.sequence), tag)], primer.name)
    row = strand_row(shown, "bottom", indent=gutter + hit.start,
                     prefix="3'- ", suffix=" -5'")
    strand = row.visual
    assert isinstance(strand, StrandVisual)
    row.visual = (strand, ArrowVisual(strand.start - 10, strand.start - 2,
                                      "left", primer.name))
    return row


@dataclass(frozen=True)
class PCRGeometry:
    template: Construct | Scene
    forward: PCRPrimer | None
    reverse: PCRPrimer | None
    availability: str
    statement: str
    forward_hit: seqprimers.Landing | None
    reverse_hit: seqprimers.Landing | None
    reverse_region: str = ""
    reverse_pool: str = ""
    symmetric_species: bool = False

    def __post_init__(self) -> None:
        if self.availability not in ("exact", "regional", "unavailable"):
            raise ValueError("PCR geometry must be exact, regional or unavailable")
        if not self.statement.strip():
            raise ValueError("PCR geometry needs a source-scoped statement")
        if isinstance(self.template, Scene) and self.availability != "unavailable":
            raise ValueError("a Scene PCR template supports unavailable geometry only")
        if self.availability == "exact":
            if (self.forward is None or self.forward_hit is None
                    or self.reverse is None or self.reverse_hit is None):
                raise ValueError("exact PCR geometry needs both primer sites")
            if self.forward_hit.strand != "bottom":
                raise ValueError("forward PCR primer must anneal to the bottom strand")
            if self.reverse_hit.strand != "top":
                raise ValueError("reverse PCR primer must anneal to the top strand")
            if self.forward_hit.end > self.reverse_hit.start:
                raise ValueError("PCR primers do not face inward around an amplicon")
        elif self.availability == "regional":
            if self.forward is None or self.forward_hit is None:
                raise ValueError("regional PCR geometry needs an exact forward site")
            if self.forward_hit.strand != "bottom":
                raise ValueError("forward PCR primer must anneal to the bottom strand")
            if self.reverse is not None or self.reverse_hit is not None:
                raise ValueError("regional geometry uses a named pool, not a false exact site")
            self.template.span(self.reverse_region)
        elif (self.reverse is not None or self.reverse_hit is not None
              or self.reverse_region or self.reverse_pool):
            raise ValueError("unavailable PCR geometry cannot claim a reverse site")

    @classmethod
    def exact(cls, template: Construct, forward: PCRPrimer, reverse: PCRPrimer,
              *, statement: str) -> "PCRGeometry":
        return cls(template, forward, reverse, "exact", statement,
                   _landing(template, forward), _landing(template, reverse))

    @classmethod
    def symmetric(cls, template: Construct, primer: PCRPrimer, *,
                  statement: str) -> "PCRGeometry":
        """Exact inward sites for one physical primer species used at both ends.

        SMART/ISPCR and similar reactions use the same handle on both cDNA ends.  The
        reagent must be represented once even though its sequence has one forward-facing
        and one reverse-facing site on the template.
        """
        return cls(template, primer, primer, "exact", statement,
                   _landing_on(template, primer, "bottom"),
                   _landing_on(template, primer, "top"),symmetric_species=True)

    @classmethod
    def regional(cls, template: Construct, forward: PCRPrimer, *, reverse_pool: str,
                 reverse_region: str, statement: str) -> "PCRGeometry":
        if not reverse_pool.strip():
            raise ValueError("regional PCR geometry needs a reverse-primer pool name")
        return cls(template, forward, None, "regional", statement,
                   _landing(template, forward), None, reverse_region, reverse_pool)

    @classmethod
    def unavailable(cls, template: Construct | Scene, *, reason: str,
                    symmetric_species: bool = False) -> "PCRGeometry":
        """Record a PCR without declared usable primer/site information."""
        return cls(template, None, None, "unavailable", reason, None, None,
                   symmetric_species=symmetric_species)

    def rows(self) -> list[Row]:
        if isinstance(self.template, Scene):
            # An unassigned primer site must not manufacture a complementary
            # strand on an explicitly single-stranded precursor.
            return [*self.template.rows(), Row(chunks=[(
                "PCR primer sites unavailable — " + self.statement, None, False)])]
        label = "PCR template  "
        gutter = len(label)
        left_pad = (max(0, self.forward_hit.free5 - self.forward_hit.start)
                    if self.forward_hit is not None else 0)
        top = strand_row(self.template, "top", indent=left_pad,
                         prefix=label + "5'- ", suffix=" -3'")
        bottom = strand_row(self.template, "bottom",
                            indent=left_pad, prefix=" " * gutter + "3'- ", suffix=" -5'")
        rows = []
        if self.forward is not None and self.forward_hit is not None:
            rows.append(_primer_row(self.template, self.forward, self.forward_hit,
                                    gutter=gutter + left_pad))
        rows.extend((top, bottom))
        if self.availability == "exact":
            assert self.reverse is not None and self.reverse_hit is not None
            rows.append(_primer_row(self.template, self.reverse, self.reverse_hit,
                                    gutter=gutter + left_pad))
            start = self.forward_hit.start - self.forward_hit.free5
            end = self.reverse_hit.end + self.reverse_hit.free5
            sequence_col = left_pad + gutter + len("5'- ")
            rows.append(Row(chunks=[(" " * (sequence_col + start)
                                      + " " * max(1, end - start), None, False)],
                            visual=SpanVisual(sequence_col + start,
                                              sequence_col + end,
                                              "resulting amplicon boundaries")))
        elif self.availability == "regional":
            start, end = self.template.span(self.reverse_region)
            base = left_pad + gutter + len("5'- ") + start
            pool = f"[{self.reverse_pool}]"
            rows.append(Row(chunks=[(" " * base + pool, "r2", False)],
                            visual=(StrandVisual(base, base + len(pool), "left"),
                                    ArrowVisual(max(0, base - 10), max(1, base - 2),
                                                "left", self.reverse_pool))))
            rows.append(Row(chunks=[(" " * base + " " * max(1, end - start)
                                      + " exact transcript coordinates not reported",
                                      None, False)],
                            visual=SpanVisual(base, gutter + len("5'- ") + end,
                                              "reverse primer pool region")))
        else:
            rows.append(Row(chunks=[("reverse primer site unavailable — " + self.statement,
                                     None, False)]))
        return rows
