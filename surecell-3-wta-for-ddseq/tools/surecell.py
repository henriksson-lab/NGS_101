"""Black-box molecular model for the commercial SureCell WTA 3' kit.

Illumina's Reference Guide #1000000021452 v01 is authoritative for the reaction flow,
but prints no oligo sequences.  Consequently every reagent-derived molecular region is
created by the inferred-placeholder helper below; no secondary reconstruction enters
the model as established sequence.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, feature


def inferred(name: str, token: str, tag: str | None = None, **kw) -> Segment:
    """Create a visibly uncertain, non-sequence token for a proprietary kit region."""
    if not token.startswith("[") or not token.endswith("]"):
        raise ValueError("proprietary region tokens must be bracketed descriptions")
    return Segment(name=name, top=token, tag=tag, placeholder=True, inferred=True,
                   note="sequence and exact length not published by Illumina", **kw)


def inferred_poly_t(name: str = "oligo-dT", display_bases: int = 12) -> Segment:
    """Draw the documented poly(dT) role without asserting its unpublished length."""
    if display_bases < 2:
        raise ValueError("poly(dT) display must contain enough bases to be legible")
    return Segment(name, "T" * display_bases, inferred=True,
                   note="poly(dT) role is documented; exact length is not published")


def bead_primer() -> Construct:
    return Construct([
        inferred("3' Barcode Mix oligo", "[proprietary barcode region]", "cbc",
                 feature=feature("cell_barcode_region", "cell_barcode", "unknown",
                                 note="internal barcode anatomy and length are proprietary")),
        inferred_poly_t(),
    ], name="SureCell bead oligo (schematic, not to scale)")


def first_strand() -> Construct:
    return Construct([
        *bead_primer().segments,
        Segment("antisense cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
    ], name="SureCell first-strand cDNA")


def selected_fragment() -> Construct:
    """Bead-end fragment selected by TPP1 and one N7xx adapter during PCR."""
    return Construct([
        inferred("TPP1-selected end", "[TPP1-selected barcode end]", "cbc",
                 feature=feature("cell_barcode_region", "cell_barcode", "unknown",
                                 note="internal barcode anatomy and length are proprietary")),
        inferred_poly_t("poly(T)"),
        Segment("3' cDNA fragment", "XXXXXXXX...XXXXXXXX", placeholder=True),
        inferred("SureCell transposome end", "[transposome-derived end]", "me"),
        inferred("N7xx adapter", "[N7xx adapter, including index]", "p7",
                 feature=feature("sample_i7", "sample_index", "whitelist",
                                 whitelist="Illumina N7xx index set")),
    ], name="SureCell final library (schematic, not to scale)")


def schematic_duplex(con: Construct) -> Scene:
    """Pair known bases and explicitly exempt only non-sequence role placeholders."""
    unpaired = [s.name + "'" for s in con if s.placeholder]
    return Scene.duplex(list(con), unpaired=unpaired)


def read_boundaries() -> tuple[tuple[str, str], ...]:
    """Return only sequencing roles stated by the vendor guide."""
    con = selected_fragment()
    names = [s.name for s in con]
    required = ["TPP1-selected end", "3' cDNA fragment", "N7xx adapter"]
    if any(x not in names for x in required):
        raise ValueError("SureCell construct lacks a required read-layout boundary")
    return (
        ("Read 1", "supplied custom Sequencing Primer SP; molecular start undisclosed"),
        ("Index 1", "the selected N7xx adapter index"),
    )
