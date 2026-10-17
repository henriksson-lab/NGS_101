"""Commercial synthetic-tag geometry for GeoMx DSP NGS readout.

The vendor publishes the region order and every length, but not the SPR, UMI,
RTS-ID, Seq Code index or primer sequences.  The final 33/34-nt read arms and
29/24-nt flow-cell arms have canonical TruSeq geometry; those canonical bases
are deliberately marked inferred so the page can compute primer landing sites
without presenting proprietary sequences as transcribed facts.
"""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, Scene, Segment, feature, revcomp


UMI = feature("geomx_umi", "umi", "random")
TARGET_ID = feature(
    "geomx_rts_id", "feature_barcode", "whitelist",
    whitelist="assay-specific GeoMx Probe Kit Configuration (PKC) file")
AOI_I5 = feature(
    "geomx_aoi_i5", "spatial_barcode", "whitelist",
    group="geomx_aoi", part="i5",
    whitelist="run-specific GeoMx SeqCodeIndices.csv")
AOI_I7 = feature(
    "geomx_aoi_i7", "spatial_barcode", "whitelist",
    group="geomx_aoi", part="i7",
    whitelist="run-specific GeoMx SeqCodeIndices.csv")


def _unknown(name: str, letter: str, length: int, tag=None, *, feature=None):
    return Segment(name, letter * length, tag, placeholder=True, feature=feature,
                   length_bp=length)


def photocleaved_tag() -> Construct:
    """The vendor-published 66-nt tag order before Seq Code PCR."""
    return Construct([
        _unknown("proprietary SPR1 PCR handle", "S", 20, "r1"),
        _unknown("14-nt UMI", "U", 14, "umi", feature=UMI),
        _unknown("12-nt RTS target identifier", "T", 12, "cbc", feature=TARGET_ID),
        _unknown("proprietary SPR2 PCR handle", "Q", 20, "r2"),
    ], name="photocleaved GeoMx tag")


def library():
    """Final vendor geometry with inferred canonical bases and validated primers."""
    lib = Construct([
        seg("P5", il.P5, "p5", inferred=True),
        _unknown("i5 AOI index", "J", 8, "cbc", feature=AOI_I5),
        seg("SPR1 / TruSeq Read 1 arm", il.TRUSEQ_READ1, "r1", inferred=True),
        _unknown("14-nt UMI", "U", 14, "umi", feature=UMI),
        _unknown("12-nt RTS target identifier", "T", 12, "cbc", feature=TARGET_ID),
        # The manual labels this completed region SPR2, 34 nt.  A canonical TruSeq
        # arm is dA plus the 33-nt Index-1 site, exactly revcomp(Read 2).
        seg("SPR2 / TruSeq Read 2 arm", revcomp(il.TRUSEQ_READ2), "r2", inferred=True),
        _unknown("i7 AOI index reverse complement", "I", 8, "cbc", feature=AOI_I7),
        seg("P7 reverse complement", il.P7_RC, "p7", inferred=True),
    ], name="GeoMx DSP Seq Code library")
    primers = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["I2"], sp.TRUSEQ["R2"])
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("invalid inferred GeoMx DSP library: " + "; ".join(problems))
    return lib, primers


def probe_rows():
    sc = Scene()
    sc.strand("probe", [
        Segment("target-bound antibody or RNA probe", "[target affinity reagent]",
                placeholder=True),
        Segment("UV-photocleavable linker", "[PC linker]", placeholder=True),
        *list(photocleaved_tag()),
    ], label="probe on fixed tissue")
    sc.labels("probe")
    return sc.rows()


def final_rows(lib: Construct):
    sc = Scene.duplex(list(lib), label="GeoMx DSP library")
    sc.labels("top")
    return sc.rows()
