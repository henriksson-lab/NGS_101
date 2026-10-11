"""Primary template-directed amplification (PTA) molecular primitives.

PTA is a phi29 strand-displacement reaction whose defining operation is
incorporation of alpha-thio-dideoxynucleotides. The exonuclease-resistant
terminated products cannot be rescued and repeatedly re-amplified, shifting
the reaction from MDA-like exponential branching toward primary-template
copying. Protocols share this operation while supplying their own known or
proprietary random-primer description and downstream library preparation.
"""
from __future__ import annotations

from chemdraw import Construct, Scene, Segment


def amplification_product(*, primer_len: int | None = None,
                          primer_note: str = "exonuclease-resistant random primer",
                          insert_len: int = 40) -> Construct:
    """One displaced PTA product, with chemically constrained termini.

    ``primer_len=None`` preserves an undisclosed primer as a role token. A
    protocol with an authoritative length can opt into a concrete N-mer.
    """
    if primer_len is not None and primer_len < 1:
        raise ValueError("a known PTA random-primer length must be positive")
    primer = (Segment("random-primer-derived 5′ end", "N" * primer_len, "tso",
                      placeholder=True, length_bp=primer_len, note=primer_note)
              if primer_len is not None else
              Segment("random-primer-derived 5′ end", "[random primer]", "tso",
                      placeholder=True, note=primer_note))
    return Construct([
        primer,
        Segment("primary-template copy", "X" * insert_len, placeholder=True),
        Segment("irreversible terminator", "ddN", "w1", placeholder=True,
                note=("alpha-thio-ddN: no 3′-OH; its phosphorothioate bond resists "
                      "phi29 proofreading exonuclease")),
    ], name="terminated PTA amplicon")


def priming_rows(*, primer_len: int | None = None):
    """Random priming and strand-displacing extension on denatured genomic DNA."""
    sc = Scene()
    sc.strand("template", [Segment("genomic template", "X" * 48, placeholder=True)],
              label="primary genomic template", rev=True)
    if primer_len is None:
        sc.strand("primer", [Segment("exonuclease-resistant random primer",
                                     "[random primer]", "tso", placeholder=True)],
                  label="random primer")
    else:
        sc.anneal("primer", [Segment("exonuclease-resistant random primer",
                                     "N" * primer_len, "tso", placeholder=True,
                                     length_bp=primer_len)],
                  to="template", pair=("exonuclease-resistant random primer",
                                        "genomic template"), label="random primer")
    sc.arrow("primer", "phi29 extends with strand displacement")
    return sc.rows()


def product_rows(*, primer_len: int | None = None, primer_note: str = ""):
    con = amplification_product(primer_len=primer_len, primer_note=primer_note)
    sc = Scene()
    sc.strand("product", list(con), label="displaced PTA product")
    sc.mark("product", "irreversible terminator",
            "alpha-thio-ddN: extension stops; no 3′-OH")
    sc.labels("product")
    return sc.rows()
