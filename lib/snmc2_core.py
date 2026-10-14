"""Reusable, source-backed snmC-seq2 random-primer and PCR endpoint."""
from __future__ import annotations

from batch_ngs import seg, truseq_library
from chemdraw import feature
from pcr import PCRGeometry, PCRPrimer


INLINE = "CGATGT"  # representative P5L_AD002 random-primer barcode
RANDOM_PRIMER = "TTCCCTACACGACGCTCTTCCGATCT" + INLINE + "N" * 9
P5_INDEX_1 = "AATGATACGGCGACCACCGAGATCTACACACGATCAGACACTCTTTCCCTACACGACGCTCT"
P7_INDEX_1 = "CAAGCAGAAGACGGCATACGAGATAGGCAATGGTGACTGGAGTTCAGACGTGTGCTCTT"
P5_INDEX_1_CODE = "ACGATCAGAC"
P7_INDEX_1_CODE = "CATTGCCTAT"


def library(insert_name: str, *, insert_nt: int, library_name: str,
            feature_id: str = "cell_inline",
            whitelist: str = "published random-primer barcode set"):
    """Build a snmC-seq2-family endpoint; only Adaptase-added bases are unknown."""
    inserts = [
        seg("6-nt inline cell barcode", INLINE, "cbc",
            feature=feature(feature_id, "cell_barcode", "whitelist", whitelist=whitelist)),
        seg("random H9 (A/C/T only)", "H" * 9, placeholder=True),
        seg(insert_name, "X" * insert_nt, placeholder=True),
        seg("Adaptase-added bases", "[Adaptase bases]", placeholder=True, inferred=True,
            note="proprietary Adaptase sequence is not printed"),
    ]
    return truseq_library(inserts, library_name, dual_index=True)


def representative_pcr_geometry(final_library):
    """Exact Supplementary Table 2 pair 1 on a compatible endpoint."""
    return PCRGeometry.exact(
        final_library,
        PCRPrimer("dual-index P5 primer 1", P5_INDEX_1, anneal_nt=25),
        PCRPrimer("dual-index P7 primer 1", P7_INDEX_1, anneal_nt=27),
        statement=("snmC-seq2 Supplementary Table 2 representative pair 1; "
                   "10-nt P5 and P7 indexes."))
