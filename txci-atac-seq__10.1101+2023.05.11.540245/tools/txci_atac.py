"""Molecular construct model for txci-ATAC-seq."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, revcomp


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


TN5_BARCODE = "GAACCGCG"       # representative A1 well, Supplementary Table S5
I7_OLIGO_INDEX = "TCGCCTTA"    # P7.S701, Supplementary Table S7
SHORT_SBS = il.TRUSEQ_READ2[-18:]


def tn5_a() -> list[Segment]:
    return [seg("s5", nx.S5, "s5"), seg("mosaic end", nx.ME, "me")]


def tn5_b() -> list[Segment]:
    return [seg("partial TruSeq Read 2", SHORT_SBS, "r2"),
            seg("Tn5 barcode", TN5_BARCODE, "cbc"),
            seg("mosaic end", nx.ME, "me")]


def tn5_bottom() -> list[Segment]:
    return [seg("mosaic-end reverse complement", nx.ME_RC, "me")]


def transposome_scene(which: str = "A") -> Scene:
    transfer = tn5_a() if which == "A" else tn5_b()
    sc = Scene()
    sc.strand("transfer", transfer, label=f"Tn5ME-{which}")
    sc.anneal("bottom", tn5_bottom(), to="transfer",
              pair=("mosaic-end reverse complement", "mosaic end"),
              label="Tn5MErev", mod5="p")
    return sc


def bead_oligo() -> list[Segment]:
    return [seg("P5", il.P5, "p5"),
            seg("GEM barcode", "N" * 16, "cbc", placeholder=True),
            seg("s5", nx.S5, "s5")]


def p7_primer() -> list[Segment]:
    return [seg("P7", il.P7, "p7"), seg("i7 oligo index", I7_OLIGO_INDEX, "cbc"),
            seg("TruSeq Read 2", il.TRUSEQ_READ2, "r2")]


SEQ_PRIMERS = (sp.NEXTERA["R1"], sp.TRUSEQ["I1"],
               sp.NEXTERA["I2"], sp.TRUSEQ["R2"])


def final_library(insert_nt: int = 32) -> Construct:
    lib = Construct([
        *bead_oligo(), seg("mosaic end", nx.ME, "me"),
        seg("accessible genomic DNA", "X" * insert_nt, placeholder=True),
        seg("opposite mosaic end", nx.ME_RC, "me"),
        seg("Tn5 barcode, read orientation", revcomp(TN5_BARCODE), "cbc"),
        seg("TruSeq Read 2 reverse complement", revcomp(il.TRUSEQ_READ2), "r2"),
        seg("i7, read orientation", revcomp(I7_OLIGO_INDEX), "cbc"),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="txci-ATAC sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid final txci-ATAC construct: " + "; ".join(problems))
    return lib


def primer_landings(lib: Construct | None = None):
    lib = lib or final_library()
    out = []
    for primer in SEQ_PRIMERS:
        hit = sp.locate(lib, primer)
        if hit is None:
            raise ValueError(f"{primer.name} does not land on {lib.name}")
        out.append((primer, hit))
    return out


def read2_layout() -> list[tuple[str, str]]:
    return [("1&ndash;8", "Tn5 barcode"), ("9&ndash;27", "mosaic end"),
            ("28&ndash;", "accessible genomic DNA")]
