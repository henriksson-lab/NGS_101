"""Small drawing helpers for multi-omic and end-mapping protocols.

The former chemistry-bearing shortcuts (contact products, capped cDNA, clicked
nascent RNA, indexed molecules and overload/unpack rows) were retired on
2026-10-19 after their last caller, meta3C, gained a native graph. Protocol
modules own their exact oligos, molecular states and reaction order.
"""
from __future__ import annotations

from chemdraw import Row, Scene, Segment


def molecule(name: str, bases: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name, bases, tag, **kw)


def analyte_split(*branches: str, source: str = "one cell") -> list[Row]:
    """A physical split: branches share provenance, not a covalent barcode."""
    if len(branches) < 2:
        raise ValueError("an analyte split needs at least two branches")
    rows = [Row(chunks=[(source, None, False)])]
    for i, branch in enumerate(branches):
        rows.append(Row(chunks=[(("├─ " if i + 1 < len(branches) else "└─ ") + branch,
                                 None, False)]))
    return rows


def native_three_prime_end(*, captured: str = "native DNA 3-prime OH") -> Scene:
    """A native DNA end joined to an end-specific adapter before bulk fragmentation."""
    parts = [molecule("genomic DNA upstream", "X" * 28, placeholder=True),
             molecule(captured, "E", placeholder=True),
             molecule("end-capture adapter", "A" * 18, "r1", placeholder=True)]
    sc = Scene(); sc.strand("captured end", parts, label="native-end capture product")
    sc.junction("captured end", captured, "end-capture adapter", "end-specific ligation")
    sc.labels("captured end")
    return sc
