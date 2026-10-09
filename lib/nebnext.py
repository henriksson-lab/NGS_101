"""Shared construct policies for current NEBNext Illumina workflows."""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Segment, feature


STANDARD_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["I2"], sp.TRUSEQ["R2"])


def indexed_library(insert: Construct, name: str, *, umi_nt: int = 0,
                    inferred_adaptors: bool = False,
                    dA_junction: bool = True) -> Construct:
    """Build a dual-index TruSeq-shaped NEBNext product and verify all four reads.

    The UMI, when present, follows the 8-base i7 index in the Index 1 read, matching
    NEB's 20-cycle i7+UMI instruction. Its bases are deliberately placeholders.
    """
    insert_parts = list(insert)
    parts = [
        Segment("P5", il.P5, "p5", inferred=inferred_adaptors),
        Segment("i5", "I" * 8, "cbc", placeholder=True, inferred=inferred_adaptors,
                feature=feature("sample_index_i5", "sample_index", "unknown")),
        Segment("Read 1 arm", il.TRUSEQ_READ1, "r1", inferred=inferred_adaptors),
        *insert_parts,
        Segment("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2", inferred=inferred_adaptors),
        Segment("i7 reverse complement", "I" * 8, "cbc", placeholder=True,
                inferred=inferred_adaptors,
                feature=feature("sample_index_i7", "sample_index", "unknown")),
    ]
    if dA_junction:
        parts.insert(3 + len(insert_parts), Segment("dA junction", "A",
                                             inferred=inferred_adaptors))
    if umi_nt:
        parts.append(Segment("UMI in i7 index read", "N" * umi_nt, "umi",
                             placeholder=True, inferred=inferred_adaptors,
                             feature=feature("umi_i7", "umi", "random")))
    parts.append(Segment("P7 reverse complement", il.P7_RC, "p7",
                         inferred=inferred_adaptors))
    lib = Construct(parts, name=name)
    errors = sp.verify(lib, STANDARD_PRIMERS, required_roles=sp.ROLES)
    if errors:
        raise ValueError("invalid NEBNext library: " + "; ".join(errors))
    return lib


def insert(name: str, nt: int = 42, **kw) -> Construct:
    return Construct([Segment(name, "X" * nt, placeholder=True, **kw)], name=name)
