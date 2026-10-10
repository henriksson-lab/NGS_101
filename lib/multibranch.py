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
                    encoding: str = "unknown"):
    """Build and primer-validate a final branch around structured identifiers."""
    insert = [*identifier_segments(cell=cell, cell_parts=cell_parts, umi=umi,
                                   guide=guide, target=target, encoding=encoding),
              seg(payload, "X" * 34, placeholder=True)]
    if architecture == "truseq":
        return truseq_library(insert, name)
    if architecture == "nextera":
        return nextera_library(insert, name)
    raise ValueError(f"unsupported Illumina architecture {architecture!r}")
