"""Molecular construct model for 10x Chromium Single Cell ATAC v1."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, complement_segments, feature

CELL_BARCODE = feature("cell_barcode", "cell_barcode", "whitelist",
                       whitelist="10x-chromium-atac-v1")
SAMPLE_INDEX = feature("sample_index_i7", "sample_index", "fixed")
FEATURES = {"10x cell barcode": CELL_BARCODE, "i7 sample index read": SAMPLE_INDEX}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


I7_INDEX_READ = "TAAGGCGA"


def bead_oligo() -> list[Segment]:
    return [seg("P5", il.P5, "p5"),
            seg("10x cell barcode", "B" * 16, "cbc", placeholder=True),
            seg("s5", nx.S5, "s5")]


def transposed_scene(insert_nt: int = 28) -> Scene:
    gap = nx.TAGMENTATION_GAP
    core = seg("accessible genomic DNA", "X" * insert_nt, placeholder=True)
    top = [seg("s5", nx.S5, "s5"), seg("left ME", nx.ME, "me"),
           seg("left 9-nt gap", "X" * gap, placeholder=True), core]
    bottom = [seg("s7", nx.S7, "s7"), seg("right ME", nx.ME, "me"),
              seg("right 9-nt gap", "X" * gap, placeholder=True),
              *complement_segments([core])]
    sc = Scene()
    sc.strand("top", top, label="DNA")
    sc.anneal("bottom", bottom, to="top", pair=("accessible genomic DNA'", "accessible genomic DNA"), label="DNA")
    sc.anneal("left ME bottom", [seg("left ME'", nx.ME_RC, "me")], to="top",
              pair=("left ME'", "left ME"), label="", mod5="p")
    sc.anneal("right ME bottom", [seg("right ME'", nx.ME_RC, "me")], to="bottom",
              pair=("right ME'", "right ME"), label="", mod5="p", above=True)
    return sc


def gap_filled_template(insert_nt: int = 28) -> Construct:
    return Construct([seg("s7", nx.S7, "s7"), seg("mosaic end", nx.ME, "me"),
                      seg("accessible genomic DNA", "X" * insert_nt, placeholder=True),
                      seg("opposite ME reverse complement", nx.ME_RC, "me"),
                      seg("s5 reverse complement", nx.S5_RC, "s5")],
                     name="gap-filled 10x ATAC fragment")


def capture_scene(insert_nt: int = 28) -> Scene:
    template = gap_filled_template(insert_nt)
    sc = Scene()
    sc.strand("template", list(template), label="template")
    sc.anneal("bead primer", bead_oligo(), to="template", pair=("s5", "s5 reverse complement"),
              label="released gel-bead oligo")
    sc.arrow("bead primer", "linear extension through the accessible fragment")
    return sc


SEQ_PRIMERS = tuple(sp.NEXTERA[k] for k in ("R1", "I1", "I2", "R2"))


def final_library(insert_nt: int = 28) -> Construct:
    lib = Construct([
        *bead_oligo(), seg("mosaic end", nx.ME, "me"),
        seg("accessible genomic DNA", "X" * insert_nt, placeholder=True),
        seg("opposite ME reverse complement", nx.ME_RC, "me"),
        seg("s7 reverse complement", nx.S7_RC, "s7"),
        seg("i7 sample index read", I7_INDEX_READ, "cbc"),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="10x Chromium Single Cell ATAC library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid 10x ATAC library: " + "; ".join(problems))
    return lib


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")
