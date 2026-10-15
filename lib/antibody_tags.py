"""Shared molecular states for oligo-tagged antibody multimodal assays."""
from __future__ import annotations

from batch_ngs import seg, truseq_library
from chemdraw import feature


def antibody_tag(*, name: str, barcode_len: int, umi_len: int,
                 capture_len: int, cell_barcode_len: int,
                 barcode_id: str, cell_id: str,
                 bead_whitelist: str, panel_whitelist: str) -> list:
    """An antibody oligo after extension has made it library-amplifiable."""
    return [
        seg("Read 1/capture handle", "H" * capture_len, "r1", placeholder=True),
        seg("cell barcode", "C" * cell_barcode_len, "cbc", placeholder=True,
            feature=feature(cell_id, "cell_barcode", "whitelist",
                            whitelist=bead_whitelist)),
        seg("antibody-tag UMI", "U" * umi_len, "umi", placeholder=True,
            feature=feature(barcode_id + "_umi", "umi", "random")),
        seg(name, "B" * barcode_len, "cbc", placeholder=True,
            feature=feature(barcode_id, "feature_barcode", "whitelist",
                            whitelist=panel_whitelist)),
    ]


def feature_library(name: str, *, barcode_id: str, cell_id: str,
                    barcode_len: int, umi_len: int, capture_len: int,
                    cell_barcode_len: int, bead_whitelist: str,
                    panel_whitelist: str):
    return truseq_library(
        antibody_tag(name="antibody feature barcode", barcode_id=barcode_id,
                     cell_id=cell_id, barcode_len=barcode_len,
                     umi_len=umi_len, capture_len=capture_len,
                     cell_barcode_len=cell_barcode_len,
                     bead_whitelist=bead_whitelist,
                     panel_whitelist=panel_whitelist), name)
