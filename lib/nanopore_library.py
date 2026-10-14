"""Shared topology and published oligos for Oxford Nanopore libraries.

The current ONT Chemistry Technical Document publishes several adapter and cDNA-kit
oligos, but not the complete motor-bearing assembly.  Keeping those two evidence states
separate here prevents pages from either hiding published sequence or inventing the
undisclosed motor/attachment structure.
"""
from __future__ import annotations

from chemdraw import Construct, Row, Scene, Segment
from protocol_graph import SequencingPlan


# Oxford Nanopore, Chemistry Technical Document CHTD_500 v1 rev AV (24 July
# 2026), section 15.  All strings are written 5' to 3'.
LA_TOP = "CCTGTACTTCGTTCAGTTACGTATTGCT"
LA_BOTTOM = "GCAATACGTAACTGAACGAAGTACAGG"
NATIVE_ADAPTER_BOTTOM = "ACGTAACTGAACGAAGTACAGG"
RAPID_ADAPTER = LA_TOP
RTP = "CTTGCCTGTCGCTCTATCTTCAGAGGAG"
CRTA = "CTTGCGGGCGGCGGACTCTCCTCTGAAGATAGAGCGACAGGCAAG"
SSPII_HANDLE = "TTTCTGTTGGTGCTGATATTGCT"
SSPII_UMI = "TTVVVVTTVVVVTTVVVVTTVVVVTTT"
SSPII_TERMINAL = "GGG"
SSPII_AS_WRITTEN = SSPII_HANDLE + SSPII_UMI + "mGmGmG"
CDNA_PRIMER_COMMON = "ATCGCCTACCGTGACAAGAAAGTTGTCGGTGT"
CDNA_PRIMER_FORWARD = CDNA_PRIMER_COMMON + "CTTTGTGACTTGCCTGTCGCTCTATCTTC"
CDNA_PRIMER_REVERSE = CDNA_PRIMER_COMMON + "TTTCTGTTGGTGCTGATATTGC"


def role(name: str, text: str, tag: str | None = None, **kw) -> Segment:
    if not (text.startswith("[") and text.endswith("]")):
        raise ValueError("an undisclosed Nanopore component must be a bracketed role token")
    return Segment(name, text, tag, placeholder=True, **kw)


def motor_library(payload: list[Segment], *, adapter: str, name: str) -> Construct:
    """Put the motor-bearing sequencing adapter before a validated payload."""
    if not payload:
        raise ValueError("a Nanopore motor adapter needs a DNA payload")
    return Construct([
        role("motor-loaded sequencing adapter", f"[motor + {adapter}]", "r2"),
        *payload,
    ], name=name)


def single_strand_rows(library: Construct, junctions=()) -> list[Row]:
    """Draw the strand in pore-entry order and mark declared covalent boundaries."""
    sc = Scene()
    sc.strand("pore-entry strand", list(library), label="pore entry")
    for left, right, label in junctions:
        names = [s.name for s in library]
        if left not in names or right not in names or names.index(right) != names.index(left) + 1:
            raise ValueError(f"{left!r} and {right!r} are not an adjacent junction")
        sc.junction("pore-entry strand", left, right, label)
    sc.labels("pore-entry strand")
    return sc.rows()


def no_primer_ending(adapter: str) -> str:
    return (f"No synthesis primer binds. The motor carried by the {adapter} docks at the "
            "nanopore and meters the attached DNA strand through the pore, adapter-bearing "
            "end first.")


def motor_sequencing_plan(library: Construct, *, limitation: str,
                          read_uses=None) -> SequencingPlan:
    """Build the primer-free sequencing endpoint shared by current ONT kits.

    A Nanopore library is sequenceable because its leading adapter carries a motor,
    not because a synthesis primer lands on the construct.  Requiring the first segment
    to be the motor-bearing adapter prevents a generic DNA product from being declared
    sequence-ready merely by supplying prose.
    """
    if not library.segments or "motor-loaded" not in library.segments[0].name:
        raise ValueError("a Nanopore sequencing endpoint must begin with the motor adapter")
    if not limitation.strip():
        raise ValueError("a proprietary motor endpoint needs a source limitation")
    return SequencingPlan(
        library=library, primers=(), acquired_roles=("Nanopore strand",),
        read_uses=read_uses or {"Nanopore strand": "include"},
        limitation=limitation)
