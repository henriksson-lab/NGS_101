"""Reusable 10x 3-prime RNA and poly(A)-mimic feature-library endpoints."""
from __future__ import annotations
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Segment, feature


SEQ_PRIMERS = (sp.NEXTERA["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])


def _shell(payload: list[Segment], name: str) -> Construct:
    lib = Construct([
        Segment("P5", il.P5, "p5"),
        Segment("10x Read 1", nx.READ1_PRIMER, "r1"),
        Segment("cell barcode", "B" * 16, "cbc", placeholder=True,
                feature=feature("cell_barcode", "cell_barcode", "whitelist",
                                whitelist="10x Chromium 3-prime whitelist")),
        Segment("UMI", "U" * 12, "umi", placeholder=True,
                feature=feature("umi", "umi", "random")),
        *payload,
        Segment("Read 2 leading A", "A", "r2"),
        Segment("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        Segment("i7 reverse complement", "I" * 8, "cbc", placeholder=True,
                feature=feature("sample_index_i7", "sample_index", "unknown")),
        Segment("P7 reverse complement", il.P7_RC, "p7"),
    ], name=name)
    problems = sp.verify(lib, SEQ_PRIMERS,
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if problems:
        raise ValueError(f"{name}: " + "; ".join(problems))
    return lib


def transcript_library(name: str = "10x 3-prime gene-expression library") -> Construct:
    return _shell([
        Segment("poly(dT)-derived tract", "T" * 20),
        Segment("cDNA insert", "X" * 40, placeholder=True),
    ], name)


def poly_a_feature_library(feature_name: str, *, feature_len: int,
                           feature_role: str = "feature_barcode",
                           feature_id: str = "feature_barcode",
                           name: str | None = None) -> Construct:
    if feature_len < 1:
        raise ValueError("a captured feature needs positive length")
    return _shell([
        Segment("poly(dT)-derived tract", "T" * 20),
        Segment(feature_name, "F" * feature_len, "cbc", placeholder=True,
                feature=feature(feature_id, feature_role, "whitelist",
                                whitelist="protocol feature set")),
        Segment("feature handle", "X" * 18, placeholder=True),
    ], name or feature_name + " library")
