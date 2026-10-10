"""Reusable bacterial single-cell indexing and library states."""
from __future__ import annotations

from batch_ngs import seg, truseq_library
from chemdraw import Row, Scene, feature


def indexed_random_rt_scene(*, method: str, barcode_len: int = 8, umi_len: int = 8) -> Scene:
    """In-cell random-primed RT carrying a plate barcode and UMI."""
    rna = [seg("bacterial RNA", "X" * 34, placeholder=True)]
    primer = [
        seg("PCR handle", "H" * 16, "r1", placeholder=True),
        seg("plate barcode", "B" * barcode_len, "cbc", placeholder=True,
            feature=feature(method + "_plate", "cell_barcode", "combinatorial",
                            group=method + "_cell", part="plate well")),
        seg("UMI", "U" * umi_len, "umi", placeholder=True,
            feature=feature(method + "_umi", "umi", "random")),
        seg("random primer", "N" * 8, placeholder=True),
    ]
    sc = Scene(); sc.strand("RNA", rna, label="fixed, permeabilized bacterium")
    sc.anneal("RT primer", primer, to="RNA", pair=("random primer", "bacterial RNA"),
              label="well-indexed random RT primer",
              unpaired=("PCR handle", "plate barcode", "UMI", "random primer"))
    sc.arrow("RT primer", "in situ reverse transcription")
    sc.labels("RT primer")
    return sc


def droplet_barcode_rows(method: str) -> list[Row]:
    return [Row(chunks=[
        ("cDNA — [plate barcode] — [UMI] — [droplet barcode]", "cbc", False)]),
        Row(chunks=[
            ("plate + droplet indexes jointly identify the bacterium", None, False)])]


def bacterial_rna_library(method: str, name: str):
    plate = feature(method + "_plate", "cell_barcode", "combinatorial",
                    group=method + "_cell", part="plate well")
    droplet = feature(method + "_droplet", "cell_barcode", "combinatorial",
                      group=method + "_cell", part="droplet")
    umi = feature(method + "_umi", "umi", "random")
    return truseq_library([
        seg("droplet barcode", "D" * 16, "cbc", placeholder=True, feature=droplet),
        seg("plate barcode", "B" * 8, "cbc", placeholder=True, feature=plate),
        seg("UMI", "U" * 8, "umi", placeholder=True, feature=umi),
        seg("bacterial cDNA", "X" * 40, placeholder=True),
    ], name)


def indexed_mda_rows() -> list[Row]:
    return [
        Row(chunks=[("one microbe per lysis droplet → MDA amplicons", None, False)]),
        Row(chunks=[("Tn5 fragments amplified DNA and transfers mosaic ends", "me", False)]),
        Row(chunks=[("barcode-bead primer extends: bead barcode — genomic insert", "cbc", False)]),
    ]
