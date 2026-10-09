"""Reusable molecular views for targeted and microbial library preparations."""
from __future__ import annotations

from chemdraw import Construct, Row, Scene, Segment, complement_segments, feature
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
    return [
        Row(chunks=[("5'-p [ligation arm]—[12-nt molecular tag]—[backbone]—[extension arm] 3'", "umi", False)]),
        Row(chunks=[("         \\________ target gap: polymerase fill ________/", None, False)]),
        Row(chunks=[("                         ligase closes probe + copied target", None, False)]),
        Row(chunks=[("linear probe and genomic DNA --exonuclease--> removed; circle survives", None, False)]),
    ]
