"""Shared molecular states for oligo-tagged antibody multimodal assays."""
from __future__ import annotations

from batch_ngs import seg, truseq_library
from chemdraw import Scene, feature


def antibody_tag(*, name: str, barcode_len: int = 15, umi_len: int = 10,
                 capture_len: int = 18, barcode_id: str = "antibody_barcode",
                 cell_id: str = "antibody_tag_cell") -> list:
    """An antibody oligo after extension has made it library-amplifiable."""
    return [
        seg("Read 1/capture handle", "H" * capture_len, "r1", placeholder=True),
        seg("cell barcode", "C" * 16, "cbc", placeholder=True,
            feature=feature(cell_id, "cell_barcode", "whitelist",
                            whitelist="single-cell bead barcode set")),
        seg("antibody-tag UMI", "U" * umi_len, "umi", placeholder=True,
            feature=feature(barcode_id + "_umi", "umi", "random")),
        seg(name, "B" * barcode_len, "cbc", placeholder=True,
            feature=feature(barcode_id, "feature_barcode", "whitelist",
                            whitelist="published antibody-tag panel")),
    ]


def bridge_extension_scene(label: str = "antibody-derived tag") -> Scene:
    """Blocked bridge oligo transfers the droplet-bead capture handle to a tag."""
    tag = [seg("antibody barcode and UMI", "B" * 25, "cbc", placeholder=True),
           seg("bridge-complementary tail", "T" * 18, "r2", placeholder=True)]
    bridge = [seg("bead-capture handle", "H" * 18, "r1", placeholder=True),
              seg("bridge annealing region", "A" * 18, "r2", placeholder=True)]
    sc = Scene(); sc.strand("tag", tag, label=label)
    sc.anneal("bridge", bridge, to="tag",
              pair=("bridge annealing region", "bridge-complementary tail"),
              label="3-prime-blocked bridge", unpaired=("bead-capture handle",))
    sc.arrow("tag", "polymerase copies the bridge handle onto the antibody tag")
    sc.labels("tag"); sc.labels("bridge")
    return sc


def feature_library(name: str, *, barcode_id: str = "antibody_barcode",
                    cell_id: str = "antibody_tag_cell"):
    return truseq_library(
        antibody_tag(name="antibody feature barcode", barcode_id=barcode_id,
                     cell_id=cell_id), name)
