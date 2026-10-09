"""Molecular construct model for sci-RNA-seq (Cao et al., 2017)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, complement_segments, feature, revcomp

UMI_FEATURE = feature("umi", "umi", "random")
RT_BARCODE = feature("cell_barcode_rt", "cell_barcode", "combinatorial",
                     group="cell_barcode", part="RT round")
I5_FEATURE = feature("sample_index_i5", "sample_index", "unknown")
I7_FEATURE = feature("sample_index_i7", "sample_index", "unknown")
FEATURES = {"UMI": UMI_FEATURE, "RT barcode": RT_BARCODE, "i5": I5_FEATURE,
            "i7": I7_FEATURE, "i7 reverse complement": I7_FEATURE}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


HANDLE = il.TRUSEQ_READ1[-18:]
RT_TEMPLATE = HANDLE + "N" * 8 + "N" * 10 + "T" * 30 + "VN"
P5_TEMPLATE = il.P5 + "N" * 10 + il.TRUSEQ_READ1
P7_TEMPLATE = il.P7 + "N" * 10 + nx.S7


def rt_primer() -> list[Segment]:
    return [seg("Read 1 handle", HANDLE, "r1"),
            seg("UMI", "U" * 8, "umi", placeholder=True),
            seg("RT barcode", "B" * 10, "cbc", placeholder=True),
            seg("dT30", "T" * 30), seg("V", "V", placeholder=True),
            seg("N", "N", placeholder=True)]


def p5_primer() -> list[Segment]:
    return [seg("P5", il.P5, "p5"), seg("i5", "I" * 10, "cbc", placeholder=True),
            seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1")]


def p7_primer() -> list[Segment]:
    return [seg("P7", il.P7, "p7"), seg("i7", "J" * 10, "cbc", placeholder=True),
            seg("s7", nx.S7, "s7")]


def rt_scene() -> Scene:
    mrna = [seg("RNA body", "X" * 28, placeholder=True), seg("poly(A)", "A" * 30)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA")
    sc.anneal("RT primer", rt_primer(), to="mRNA", pair=("dT30", "poly(A)"),
              label="barcoded RT primer", unpaired=("V", "N"))
    sc.arrow("RT primer", "reverse transcriptase")
    return sc


def first_strand() -> Construct:
    return Construct([*rt_primer(), seg("antisense cDNA", "X" * 34, placeholder=True)],
                     name="sci-RNA-seq first strand")


def second_strand_scene() -> Scene:
    first = first_strand()
    sc = Scene.duplex(list(first), label="cDNA")
    sc.labels("top")
    return sc


def selected_tagmented_strand(insert_nt: int = 34) -> Construct:
    """RT-primer end through the nearest s7 insertion after gap fill."""
    return Construct([
        seg("Read 1 handle", HANDLE, "r1"),
        seg("UMI", "U" * 8, "umi", placeholder=True),
        seg("RT barcode", "B" * 10, "cbc", placeholder=True),
        seg("poly(T)", "T" * 30), seg("V", "V", placeholder=True),
        seg("antisense cDNA", "X" * insert_nt, placeholder=True),
        seg("mosaic end reverse complement", nx.ME_RC, "me"),
        seg("s7 reverse complement", nx.S7_RC, "s7"),
    ], name="selected sci-RNA-seq tagmented strand")


SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.NEXTERA["I1"],
               sp.TRUSEQ["I2"], sp.NEXTERA["R2"])


def final_library(insert_nt: int = 34) -> Construct:
    selected = selected_tagmented_strand(insert_nt)
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("i5", "I" * 10, "cbc", placeholder=True),
        seg("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
        *list(selected)[1:],
        seg("i7 reverse complement", "J" * 10, "cbc", placeholder=True),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="sci-RNA-seq sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid final sci-RNA-seq construct: " + "; ".join(problems))
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


def read1_layout() -> list[tuple[str, str]]:
    return [(x.cycles, x.feature.label) for x in
            sp.feature_spans(final_library(), SEQ_PRIMERS, {"Read 1": 18})]
