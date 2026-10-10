"""Reusable molecular states for non-standard RNA amplification workflows.

These helpers describe transformations of RNA itself.  They deliberately stop before
platform adapters are installed; protocol modules remain responsible for the exact final
library and sequencing-primer declarations.
"""
from __future__ import annotations

from chemdraw import Construct, Row, Scene, Segment, revcomp


def seg(name: str, bases: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name, bases, tag, **kw)


def poly_a_capture(primer: list[Segment], *, annealing_segment: str,
                   transcript_width: int = 30, poly_a_length: int | None = None,
                   primer_label: str = "capture primer") -> Scene:
    """Anneal a compound 3-prime capture primer to an RNA poly(A) tail.

    The primer's annealing segment determines the tail length, preventing a page from
    drawing a partial overlap whose geometry was maintained separately by hand.
    """
    anneal = next((s for s in primer if s.name == annealing_segment), None)
    if anneal is None:
        raise ValueError(f"capture primer has no segment {annealing_segment!r}")
    n = poly_a_length or len(anneal)
    if n != len(anneal):
        raise ValueError("the drawn poly(A) tract must equal the annealing segment length")
    rna = [seg("transcript body", "X" * transcript_width, placeholder=True),
           seg("poly(A)", "A" * n)]
    sc = Scene(); sc.strand("mRNA", rna, label="mRNA")
    sc.anneal("primer", primer, to="mRNA", pair=(annealing_segment, "poly(A)"),
              label=primer_label,
              unpaired=tuple(s.name for s in primer if s.name != annealing_segment))
    sc.labels("primer")
    return sc


def terminal_tail(template: Construct, base: str, count: int, *, enzyme: str,
                  name: str | None = None) -> tuple[Construct, Scene]:
    """Append a homopolymer tail and draw it, deriving the product from the input."""
    if base.upper() not in "ACGT" or len(base) != 1 or count < 1:
        raise ValueError("terminal_tail needs one DNA base and a positive length")
    tail = seg(name or f"poly(d{base.upper()}) tail", base.upper() * count, "tso")
    out = Construct([*template.segments, tail], name=f"{template.name} + d{base.upper()} tail")
    sc = Scene(); sc.strand("product", list(out), label="tailed single strand")
    sc.mark("product", tail.name, f"{enzyme} adds d{base.upper()}")
    sc.labels("product")
    return out, sc


def direct_ssrna_ivt(template: list[Segment], promoter: Segment) -> Scene:
    """Draw a duplex T7 promoter attached to an otherwise single-stranded RNA template."""
    sc = Scene(); sc.strand("RNA template", [promoter, *template], label="original ssRNA")
    # ``Scene.anneal`` accepts every strand 5′→3′ and reverses its drawing, so the
    # partner must be the reverse complement, not the left-to-right plain complement.
    partner = Segment("T7 promoter partner", revcomp(promoter.top))
    sc.anneal("promoter partner", [partner], to="RNA template",
              pair=(partner.name, promoter.name), label="short promoter strand")
    sc.mark("RNA template", promoter.name, "duplex T7 promoter")
    sc.arrow("RNA template", "T7 RNA polymerase transcribes along the original ssRNA")
    return sc


def linear_amplification_rows(*, source: str, product: str) -> list[Row]:
    """Proportional-font explanation beneath a molecular IVT state."""
    return [Row(chunks=[(f"one {source} → many antisense {product} copies", None, False)])]
