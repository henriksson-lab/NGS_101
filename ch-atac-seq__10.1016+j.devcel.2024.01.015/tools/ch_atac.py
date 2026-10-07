"""Conservative molecular model for CH-ATAC-seq.

The defining paper's oligo table was unavailable.  Every ATAC construct is therefore
marked inferred; constants transcribed from the related CH-RNA-seq paper are retained as
source boundaries rather than silently promoted to CH-ATAC evidence.
"""
from __future__ import annotations

import nextera as nx
import rt
import seqprimers as sp
from chemdraw import Construct, Segment, revcomp


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


MGI_P5 = "GAACGACATGGCTACGATCCGACTT"
MGI_P7 = "TGTGAGCCAAGGAGTTGTTGTCTTC"
HY_BARCODE = "TTCTCGCATG"       # representative CH-RNA HY barcode
TN5_BARCODE = "N" * 10          # CH-ATAC table unavailable
I7_OLIGO_INDEX = "TAGGTCCGAT"   # representative CH-RNA MGI P7 index


def inferred(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return seg(name, top, tag, inferred=True, **kw)


def indexed_tn5_transfer() -> list[Segment]:
    return [inferred("SMART handle", rt.SMART_HANDLE, "tso"),
            inferred("Tn5 barcode", TN5_BARCODE, "cbc", placeholder=True),
            inferred("mosaic end", nx.ME, "me")]


def primer_c() -> list[Segment]:
    return [inferred("s7", nx.S7, "s7"), inferred("mosaic end", nx.ME, "me")]


def hy_head() -> list[Segment]:
    return [inferred("s5", nx.S5, "s5"), inferred("mosaic end", nx.ME, "me")]


def barcoded_hy() -> list[Segment]:
    return [inferred("SMART-handle reverse complement", revcomp(rt.SMART_HANDLE), "tso"),
            inferred("HY barcode", HY_BARCODE, "cbc"),
            inferred("HY-head reverse complement", revcomp(nx.ADAPTOR_S5), "s5")]


SEQ_PRIMERS = (sp.NEXTERA["R1"], sp.NEXTERA["I1"], sp.NEXTERA["R2"])


def final_library(insert_nt: int = 32) -> Construct:
    lib = Construct([
        inferred("MGI P5", MGI_P5, "p5"), inferred("s5", nx.S5, "s5"),
        inferred("mosaic end", nx.ME, "me"),
        inferred("HY barcode, read orientation", revcomp(HY_BARCODE), "cbc"),
        inferred("SMART handle", rt.SMART_HANDLE, "tso"),
        inferred("Tn5 barcode", TN5_BARCODE, "cbc", placeholder=True),
        inferred("mosaic end after barcodes", nx.ME, "me"),
        seg("accessible genomic DNA", "X" * insert_nt, placeholder=True),
        inferred("opposite mosaic end", nx.ME_RC, "me"),
        inferred("s7 reverse complement", nx.S7_RC, "s7"),
        inferred("i7, read orientation", revcomp(I7_OLIGO_INDEX), "cbc"),
        inferred("MGI P7 reverse complement", revcomp(MGI_P7), "p7"),
    ], name="inferred CH-ATAC sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS,
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if problems:
        raise ValueError("invalid CH-ATAC reconstruction: " + "; ".join(problems))
    return lib


def primer_landings(lib: Construct | None = None):
    lib = lib or final_library()
    return [(p, sp.locate(lib, p)) for p in SEQ_PRIMERS]


def read1_layout() -> list[tuple[str, str]]:
    return [("1&ndash;10", "HY barcode"), ("11&ndash;33", "SMART handle — dark"),
            ("34&ndash;43", "Tn5 barcode"), ("44&ndash;62", "mosaic end — dark"),
            ("63&ndash;", "accessible genomic DNA")]
