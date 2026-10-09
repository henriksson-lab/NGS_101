"""Molecular model for CH-ATAC-seq from Supplementary Table 1."""
from __future__ import annotations

import nextera as nx
import rt
import seqprimers as sp
from chemdraw import Construct, Segment, feature, revcomp


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


MGI_P5 = "GAACGACATGGCTACGATCCGACTT"
MGI_P7 = "TGTGAGCCAAGGAGTTGTTGTCTTC"
HY_BARCODE = "TTCTCGCATG"       # representative CH-ATAC Barcoded_HY_oligo_1
TN5_BARCODE = "TCCTACCAGT"      # representative Tn5_barcode_primer_1
I7_OLIGO_INDEX = "TAGGTCCGAT"   # representative MGI_P7_index_1


def indexed_tn5_transfer() -> list[Segment]:
    return [seg("SMART handle", rt.SMART_HANDLE, "tso"),
            seg("Tn5 barcode", TN5_BARCODE, "cbc"),
            seg("mosaic end", nx.ME, "me")]


def primer_c() -> list[Segment]:
    return [seg("s7", nx.S7, "s7"), seg("mosaic end", nx.ME, "me")]


def hy_head() -> list[Segment]:
    return [seg("s5", nx.S5, "s5"), seg("mosaic end", nx.ME, "me")]


def barcoded_hy() -> list[Segment]:
    return [seg("SMART-handle reverse complement", revcomp(rt.SMART_HANDLE), "tso"),
            seg("HY barcode", HY_BARCODE, "cbc"),
            seg("HY-head reverse complement", revcomp(nx.ADAPTOR_S5), "s5")]


SEQ_PRIMERS = (sp.NEXTERA["R1"], sp.NEXTERA["I1"], sp.NEXTERA["R2"])


def final_library(insert_nt: int = 32) -> Construct:
    lib = Construct([
        seg("MGI P5", MGI_P5, "p5"), seg("s5", nx.S5, "s5"),
        seg("mosaic end", nx.ME, "me"),
        seg("HY barcode, read orientation", revcomp(HY_BARCODE), "cbc",
            feature=feature("cell_hy", "cell_barcode", "combinatorial", group="cell_id", part="HY")),
        seg("SMART handle", rt.SMART_HANDLE, "tso"),
        seg("Tn5 barcode", TN5_BARCODE, "cbc",
            feature=feature("cell_tn5", "cell_barcode", "combinatorial", group="cell_id", part="Tn5")),
        seg("mosaic end after barcodes", nx.ME, "me"),
        seg("accessible genomic DNA", "X" * insert_nt, placeholder=True),
        seg("opposite mosaic end", nx.ME_RC, "me"),
        seg("s7 reverse complement", nx.S7_RC, "s7"),
        seg("i7, read orientation", revcomp(I7_OLIGO_INDEX), "cbc",
            feature=feature("sample_i7", "sample_index", "whitelist", whitelist="published MGI P7 index set")),
        seg("MGI P7 reverse complement", revcomp(MGI_P7), "p7"),
    ], name="CH-ATAC sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS,
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if problems:
        raise ValueError("invalid CH-ATAC library: " + "; ".join(problems))
    return lib


def primer_landings(lib: Construct | None = None):
    lib = lib or final_library()
    return [(p, sp.locate(lib, p)) for p in SEQ_PRIMERS]


def read1_layout() -> list[tuple[str, str]]:
    return [("1&ndash;10", "HY barcode"), ("11&ndash;33", "SMART handle — dark"),
            ("34&ndash;43", "Tn5 barcode"), ("44&ndash;62", "mosaic end — dark"),
            ("63&ndash;", "accessible genomic DNA")]
