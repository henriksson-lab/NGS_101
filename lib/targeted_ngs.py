"""Bounded library construction for amplicons with two indexing PCR stages."""
from batch_ngs import nextera_library, seg


def two_stage_amplicon(target_name: str, target_length: int = 34):
    """Nextera-tailed locus PCR followed by dual-index PCR."""
    insert = [seg(target_name, "N" * target_length, placeholder=True)]
    return nextera_library(insert, f"dual-indexed {target_name} amplicon")
