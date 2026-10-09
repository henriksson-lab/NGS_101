"""Reusable construct models for Zhang-lineage pooled CRISPR screens.

The protocol pages choose a vector and a readout independently.  This module keeps that
separation structural: vector pages cannot silently turn a one-PCR readout into a nested
readout, and both readout builders return finished libraries on which ``seqprimers`` can
locate the primers actually used by the run.
"""
from __future__ import annotations

import crispr
import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp


def seg(name: str, seq: str, tag: str | None = None, *, placeholder: bool = False,
        **kw) -> Segment:
    return Segment(name, seq, tag, placeholder=placeholder, **kw)


SPACER = "N" * crispr.SPACER_LEN


def guide_oligos() -> tuple[str, str]:
    """BsmBI-compatible oligos for a representative non-G-starting spacer."""
    return crispr.guide_oligos(SPACER)


def cloned_guide_cassette() -> Construct:
    """The sequence restored by guide ligation, with the two ligation seams named."""
    return Construct([
        seg("hU6 promoter", crispr.U6_UPSTREAM + crispr.U6_3PRIME, "s5"),
        seg("Pol III +1 G", crispr.U6_PLUS1, "me"),
        seg("20-nt guide", SPACER, "cbc", placeholder=True,
            feature=feature("guide", "guide_barcode", "unknown")),
        seg("sgRNA scaffold", crispr.SCAFFOLD_V1, "r2"),
    ], name="cloned guide cassette")


def cloning_scene() -> Scene:
    """Duplex cassette with both BsmBI-directed ligation junctions marked."""
    con = cloned_guide_cassette()
    sc = Scene.duplex(list(con), label="cloned cassette")
    sc.junction("top", "hU6 promoter", "Pol III +1 G", "ligation junction 1")
    sc.junction("top", "20-nt guide", "sgRNA scaffold", "ligation junction 2")
    return sc


def vector_map(version: str) -> tuple[tuple[str, str], ...]:
    """Ordered transfer-genome features; cPPT placement is the v1/v2 distinction."""
    common_left = (("5' LTR", "r3"), ("packaging / RRE", "s5"))
    guide = (("hU6-guide-scaffold", "cbc"),)
    expression = (("EFS-Cas9-P2A-Puro", "me"), ("WPRE", "t7"), ("3' LTR", "r2"))
    if version == "v1":
        return common_left + guide + (("cPPT/CTS", "w1"),) + expression
    if version == "v2":
        return common_left + (("cPPT/CTS", "w1"),) + guide + expression
    raise ValueError("version must be 'v1' or 'v2'")


# Joung et al. 2017: a single PCR, with a reverse primer in the invariant scaffold.
JOUNG_STAGGER = "TAAGTAGAG"
JOUNG_U6_SITE = crispr.U6_UPSTREAM[1:] + crispr.U6_3PRIME
JOUNG_SCAFFOLD_SITE_TOP = crispr.SCAFFOLD_V1[50:73]
JOUNG_I7 = "TCGCCTTG"


def one_pcr_library() -> Construct:
    """Finished ~260-bp Joung one-PCR library for a non-G-starting guide."""
    con = Construct([
        seg("P5", il.P5, "p5"),
        seg("Read 1 site", il.TRUSEQ_P5_FULL[len(il.P5):], "s5"),
        seg("stagger", JOUNG_STAGGER, "me"),
        seg("hU6 primer site", JOUNG_U6_SITE, "s5"),
        seg("Pol III +1 G", crispr.U6_PLUS1, "me"),
        seg("20-nt guide", SPACER, "cbc", placeholder=True,
            feature=feature("guide", "guide_barcode", "unknown")),
        seg("scaffold before reverse site", crispr.SCAFFOLD_V1[:50], "r2"),
        seg("reverse-primer site", JOUNG_SCAFFOLD_SITE_TOP, "r2"),
        seg("Read 2 site", revcomp(il.TRUSEQ_READ2), "t7"),
        seg("i7", revcomp(JOUNG_I7), "cbc",
            feature=feature("sample_index_i7", "sample_index", "fixed")),
        seg("P7", revcomp(il.P7), "p7"),
    ], name="Joung one-PCR library")
    problems = sp.verify(con, one_pcr_sequencing_primers(),
                         required_roles=("Read 1", "Index 1 (i7)"))
    if problems:
        raise ValueError("invalid one-PCR library: " + "; ".join(problems))
    return con


def one_pcr_sequencing_primers():
    return [sp.TRUSEQ["R1"], sp.TRUSEQ["I1"]]


# Moffat Lab TKOv3 two-step readout (2019 protocol, LCV2::TKOv3 branch).
MOFFAT_PCR1_F = "GAGGGCCTATTTCCCATGATTC"
MOFFAT_PCR1_R = "GTTGCGAAAAAGAACGTTCACGG"
MOFFAT_PCR2_SITE_F = "TTGTGGAAAGGACGAAACACCG"
MOFFAT_PCR2_SITE_R = "ACTTGCTATTTCTAGCTCTAAAAC"
MOFFAT_I5 = il.NEBNEXT_I5_SET1["i501"]
MOFFAT_I7 = il.NEBNEXT_I7_SET1["i701"]


def two_step_intermediate() -> Construct:
    """PCR1 enrichment product; distal genomic/vector sequence is intentionally schematic."""
    return Construct([
        seg("PCR1 forward site", MOFFAT_PCR1_F, "s5"),
        seg("left vector / genomic flank", "X" * 12, "w1", placeholder=True),
        seg("hU6", crispr.U6_3PRIME + crispr.U6_PLUS1, "s5"),
        seg("20-nt guide", SPACER, "cbc", placeholder=True,
            feature=feature("guide", "guide_barcode", "unknown")),
        seg("scaffold", crispr.SCAFFOLD_V1[:24], "r2"),
        seg("right vector / genomic flank", "X" * 12, "w1", placeholder=True),
        seg("PCR1 reverse site", revcomp(MOFFAT_PCR1_R), "r2"),
    ], name="PCR1 enriched guide locus (~600 bp)")


def two_step_library() -> Construct:
    """Finished dual-index LCV2::TKOv3 PCR2 library from the printed primer sequences."""
    con = Construct([
        seg("P5", il.P5, "p5"),
        seg("i5", MOFFAT_I5, "cbc",
            feature=feature("sample_index_i5", "sample_index", "fixed")),
        seg("Read 1 site", il.TRUSEQ_READ1, "s5"),
        seg("PCR2 hU6 site", MOFFAT_PCR2_SITE_F, "s5"),
        seg("20-nt guide", SPACER, "cbc", placeholder=True,
            feature=feature("guide", "guide_barcode", "unknown")),
        seg("PCR2 scaffold site", revcomp(MOFFAT_PCR2_SITE_R), "r2"),
        seg("Read 2 site", revcomp(il.TRUSEQ_READ2), "t7"),
        seg("i7", revcomp(MOFFAT_I7), "cbc",
            feature=feature("sample_index_i7", "sample_index", "fixed")),
        seg("P7", revcomp(il.P7), "p7"),
    ], name="Moffat two-step PCR library")
    problems = sp.verify(con, two_step_sequencing_primers(),
                         required_roles=("Read 1", "Index 1 (i7)", "Index 2 (i5)"))
    if problems:
        raise ValueError("invalid two-step library: " + "; ".join(problems))
    return con


def two_step_sequencing_primers():
    return [sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["I2"]]
