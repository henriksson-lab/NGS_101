"""Molecular model for VASA-seq as defined by Salmen et al. (2022).

The defining paper prints the shared RA3/RTP pair and VASA-drop index-primer
structures.  Capture oligos and the VASA-plate PCR arms are cited to earlier
protocols but are not printed, so those regions remain inferred placeholders.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from chemdraw import Construct, Scene, Segment, revcomp
from illumina import P5, P7
import nextera as nx


RA3 = "TGGAATTCTCGGGTGCCAAGG"
RTP = "GCCTTGGCACCCGAGAATTCCA"
DROP_COMMON = "GGGTGTCGGGTGCAG"
INDEX_LEN = 8
UFI_LEN = 6
DROP_BARCODE_LEN = 8
PLATE_BARCODE_LEN = 8


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def inferred(name: str, token: str, tag: str | None = None) -> Segment:
    """Create an unpublished molecular region with mandatory inferred styling."""
    if not (token.startswith("[") and token.endswith("]")):
        raise ValueError("inferred molecular regions must be bracketed role tokens")
    return seg(name, token, tag, placeholder=True, inferred=True)


def poly_a(name: str = "added poly(A)", display_bases: int = 12) -> Segment:
    """Display the published poly(A)-tailing reaction without asserting tail length."""
    if display_bases < 2:
        raise ValueError("poly(A) display must be legible")
    return seg(name, "A" * display_bases, inferred=True)


def repaired_fragment() -> Construct:
    return Construct([
        seg("RNA fragment", "XXXXXXXX...XXXXXXXX", placeholder=True),
        poly_a(),
    ], name="fragmented, repaired and poly(A)-tailed cellular RNA")


def capture_primer(fmt: str) -> Construct:
    if fmt not in {"plate", "drop"}:
        raise ValueError("VASA format must be 'plate' or 'drop'")
    return Construct([
        inferred(f"{fmt} barcode / T7 region", f"[{fmt} capture region]", "cbc"),
        seg("UFI", "N" * UFI_LEN, "umi", placeholder=True),
        seg("cell barcode", "N" * (PLATE_BARCODE_LEN if fmt == "plate" else DROP_BARCODE_LEN),
            "cbc", placeholder=True),
        seg("oligo-dT", "T" * 12, inferred=True),
    ], name=f"VASA-{fmt} capture oligo (schematic)")


def ivt_template(fmt: str) -> Construct:
    p = capture_primer(fmt)
    return Construct([
        *p.segments,
        seg("fragment cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
    ], name=f"VASA-{fmt} IVT template (schematic)")


def schematic_duplex(con: Construct) -> Scene:
    """Pair every sequence-like column; exempt only bracketed role placeholders."""
    unpaired = [s.name + "'" for s in con if s.placeholder and s.top.startswith("[")]
    return Scene.duplex(list(con), unpaired=unpaired)


def ligated_arna() -> Construct:
    return Construct([
        seg("rRNA-depleted aRNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
        seg("RA3", RA3, "r2"),
    ], name="rRNA-depleted aRNA after RA3 ligation")


def rtp_segments() -> list[Segment]:
    """RTP assembled so its published 21-nt annealing core is guaranteed from RA3."""
    core = revcomp(RA3)
    if RTP != "G" + core:
        raise ValueError("published RTP no longer equals 5' G plus reverse-complement RA3")
    return [seg("RTP 5' G", RTP[0]), seg("RTP core", core, "r2")]


def drop_i5_primer(index: str) -> Construct:
    if len(index) != INDEX_LEN:
        raise ValueError("VASA-drop i5 index must be 8 nt")
    return Construct([seg("P5", P5, "p5"), seg("i5", index, "cbc"),
                      seg("s5", nx.S5, "s5")], name="VD_ILMN8_i5")


def drop_i7_primer(index: str) -> Construct:
    if len(index) != INDEX_LEN:
        raise ValueError("VASA-drop i7 index must be 8 nt")
    return Construct([seg("P7", P7, "p7"), seg("i7", index, "cbc"),
                      seg("15-nt common sequence", DROP_COMMON)], name="VD_ILMN8_i7")


def plate_library() -> Construct:
    """Final VASA-plate molecule at the resolution supported by the defining paper."""
    con = Construct([
        inferred("plate P5-side PCR arm", "[plate P5-side arm]", "p5"),
        seg("UFI", "N" * UFI_LEN, "umi", placeholder=True),
        seg("cell barcode", "N" * PLATE_BARCODE_LEN, "cbc", placeholder=True),
        seg("poly(T) read-through", "T" * 12),
        seg("fragment cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
        seg("RA3", RA3, "r2"),
        inferred("plate P7-side PCR arm", "[plate P7-side arm]", "p7"),
    ], name="VASA-plate final library")
    _require_order(con, ("UFI", "cell barcode", "poly(T) read-through", "fragment cDNA", "RA3"))
    return con


def drop_library() -> Construct:
    """Final VASA-drop molecule using paper-supported read boundaries.

    The disputed ligation-adapter/bead-linker sequence is deliberately represented
    only by a role placeholder.
    """
    con = Construct([
        *drop_i5_primer("N" * INDEX_LEN).segments,
        seg("fragment cDNA", "XXXXXXXX...XXXXXXXX", placeholder=True),
        poly_a("poly(A) copy"),
        seg("UFI'", "N" * UFI_LEN, "umi", placeholder=True),
        seg("barcode 2'", "N" * DROP_BARCODE_LEN, "cbc", placeholder=True),
        inferred("bead linker", "[proprietary bead linker]"),
        seg("barcode 1'", "N" * DROP_BARCODE_LEN, "cbc", placeholder=True),
        seg("common sequence'", revcomp(DROP_COMMON)),
        seg("i7'", "N" * INDEX_LEN, "cbc", placeholder=True),
        seg("P7'", revcomp(P7), "p7"),
    ], name="VASA-drop final library")
    _require_order(con, ("fragment cDNA", "UFI'", "barcode 2'", "bead linker",
                         "barcode 1'", "common sequence'", "i7'"))
    return con


def _require_order(con: Construct, ordered: tuple[str, ...]) -> None:
    names = [s.name for s in con]
    try:
        positions = [names.index(name) for name in ordered]
    except ValueError as exc:
        raise ValueError(f"{con.name} lacks a read-layout segment") from exc
    if positions != sorted(positions):
        raise ValueError(f"{con.name} read-layout segments are out of order")


def read_layouts() -> tuple[tuple[str, str, str], ...]:
    """Read boundaries reported by the defining paper, checked against the models."""
    plate, drop = plate_library(), drop_library()
    _require_order(plate, ("UFI", "cell barcode", "poly(T) read-through"))
    _require_order(drop, ("fragment cDNA", "UFI'", "barcode 2'", "barcode 1'",
                          "common sequence'", "i7'"))
    return (
        ("VASA-plate", "Read 1 · 26 cycles", "UFI 6 + cell barcode 8 + poly(T)"),
        ("VASA-plate", "Read 2 · 135 cycles", "fragment cDNA"),
        ("VASA-drop", "Read 1 · 247 cycles", "fragment cDNA"),
        ("VASA-drop", "Index 1 · 31 cycles", "barcode 1 (8) + common sequence (15) + i7 (8)"),
        ("VASA-drop", "Index 2 · 8 cycles", "i5"),
        ("VASA-drop", "Read 2 · 14 cycles", "barcode 2 (8) + UFI (6)"),
    )
