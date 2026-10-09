"""Reusable molecular views for targeted and microbial library preparations."""
from __future__ import annotations

from chemdraw import (Construct, Row, Scene, Segment, circle_rows, complement_segments,
                      feature)
from batch_ngs import nextera_library, seg, truseq_library


def hybrid_capture_scene() -> Scene:
    """An indexed library hybridized to one complementary biotinylated RNA bait."""
    library = [seg("left adapter", "X" * 12, placeholder=True),
               seg("off-target flank", "N" * 8, placeholder=True),
               seg("target interval", "ACGTTGCAACGTTGCAACGT"),
               seg("right flank / adapter", "Y" * 16, placeholder=True)]
    bait = [seg("biotinylated RNA bait", "ACGTTGCAACGTTGCAACGT")]
    sc = Scene()
    sc.strand("library", library, label="indexed library")
    sc.anneal("bait", complement_segments(bait), to="library",
              pair=("biotinylated RNA bait'", "target interval"),
              label="5'-biotin RNA bait")
    sc.mark("bait", "biotinylated RNA bait'", "streptavidin capture")
    return sc


def two_stage_amplicon(target_name: str, target_length: int = 34):
    """Nextera-tailed locus PCR followed by dual-index PCR."""
    insert = [seg(target_name, "N" * target_length, placeholder=True)]
    return nextera_library(insert, f"dual-indexed {target_name} amplicon")


def shotgun_library():
    return nextera_library([seg("metagenomic DNA fragment", "N" * 36,
                                placeholder=True)],
                           "Illumina DNA Prep shotgun library")


def sureselect_library():
    """HS2 MBC form: five inline molecular-barcode bases at each insert end."""
    return truseq_library([
        seg("left molecular barcode", "U" * 5, "umi", placeholder=True,
            feature=feature("molecular_barcode_left", "umi", "random",
                            group="molecular_barcode", part="left")),
        seg("captured target insert", "N" * 34, placeholder=True),
        seg("right molecular barcode", "V" * 5, "umi", placeholder=True,
            feature=feature("molecular_barcode_right", "umi", "random",
                            group="molecular_barcode", part="right")),
    ], "SureSelect XT HS2 captured library")


def smmip_rows() -> list[Row]:
    """Closed smMIP capture product after gap fill and ligation.

    The ordered segments make the circular junction unambiguous: polymerase extends the
    extension arm across the captured gap, then ligase joins that copy to the ligation
    arm.  ``circle_rows`` supplies strand geometry and segment-derived hover text instead
    of leaving the molecular structure embedded in punctuation.
    """
    circle = Construct([
        Segment("ligation targeting arm", "X" * 18, placeholder=True),
        Segment("12-nt molecular tag", "U" * 12, "umi", placeholder=True,
                feature=feature("molecular_tag", "umi", "random")),
        Segment("invariant probe backbone", "N" * 24, placeholder=True),
        Segment("extension targeting arm", "Y" * 18, placeholder=True),
        Segment("copied target gap", "N" * 28, placeholder=True),
    ], name="closed smMIP capture product")
    return circle_rows(circle, "gap fill + ligation")
