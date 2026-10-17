"""Shared final-library constructors for protocols with several molecular branches."""
from __future__ import annotations

from batch_ngs import nextera_library, seg, truseq_library
from chemdraw import feature


def identifier_segments(*, cell: int = 0,
                        cell_parts: tuple[tuple[str, int], ...] = (),
                        umi: int = 0, guide: int = 0,
                        target: int = 0, encoding: str = "unknown"):
    out = []
    if cell and cell_parts:
        raise ValueError("choose either one cell barcode or named cell-barcode parts")
    if cell:
        out.append(seg("cell barcode", "B" * cell, "cbc", placeholder=True,
                       feature=feature("cell_barcode", "cell_barcode", encoding)))
    for part, length in cell_parts:
        key = part.lower().replace(" ", "_").replace("-", "_")
        out.append(seg(part, "B" * length, "cbc", placeholder=True,
                       feature=feature(f"cell_barcode_{key}", "cell_barcode",
                                       "combinatorial", group="cell_barcode",
                                       part=part)))
    if umi:
        out.append(seg("UMI", "U" * umi, "umi", placeholder=True,
                       feature=feature("umi", "umi", "random")))
    if guide:
        out.append(seg("guide identity", "G" * guide, "cbc", placeholder=True,
                       feature=feature("guide_barcode", "guide_barcode", encoding)))
    if target:
        out.append(seg("target identity", "T" * target, "cbc", placeholder=True,
                       feature=feature("target_barcode", "feature_barcode", encoding)))
    return out


def illumina_branch(name: str, payload: str, *, architecture: str = "truseq",
                    cell: int = 0,
                    cell_parts: tuple[tuple[str, int], ...] = (),
                    umi: int = 0, guide: int = 0, target: int = 0,
                    encoding: str = "unknown", inferred_adapters: bool = False,
                    identifier_side: str = "read1"):
    """Build and primer-validate a final branch around structured identifiers.

    The identifiers are one contiguous block on one side of the payload, and
    ``identifier_side`` says which: ``"read1"`` puts them before the payload (next
    to the Read 1 adapter), ``"read2"`` after it (next to the Read 2 adapter).  A
    protocol whose identifiers are split between reads, or sit in index reads, does
    not fit this constructor and must assemble its own library.  Only the two
    canonical Illumina shells are offered; a non-Illumina platform is out of scope.
    """
    ids = identifier_segments(cell=cell, cell_parts=cell_parts, umi=umi,
                              guide=guide, target=target, encoding=encoding)
    body = seg(payload, "X" * 34, placeholder=True)
    if identifier_side == "read1":
        insert = [*ids, body]
    elif identifier_side == "read2":
        insert = [body, *ids]
    else:
        raise ValueError(f"identifier_side must be 'read1' or 'read2', not {identifier_side!r}")
    if architecture == "truseq":
        return truseq_library(insert, name, inferred_adapters=inferred_adapters)
    if architecture == "nextera":
        return nextera_library(insert, name, inferred_adapters=inferred_adapters)
    raise ValueError(f"unsupported Illumina architecture {architecture!r}")
