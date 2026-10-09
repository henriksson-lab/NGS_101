"""Molecular construct model for Quartz-Seq2 (Sasagawa et al., 2018)."""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp

CELL_BARCODE = feature("cell_barcode", "cell_barcode", "whitelist",
                       whitelist="quartz-seq2-v3.2-rt-primer-set")
UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "fixed")
FEATURES = {"cell barcode": CELL_BARCODE, "UMI": UMI_FEATURE, "i7": I7_FEATURE}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


M = "TATAGAATTCGCGGCCGCTCGCGAT"
TAGGING_PRIMER = M + "T" * 23
GM_PRIMER = "G" + M
V32_BARCODE = "ATCAATCCTATCTGC"       # eMDRT0001
RT_PRIMER_1536_FIRST = M + "AC" + V32_BARCODE + "N" * 8 + "T" * 24
RYSHAPE_P5 = "GATCGGAAGAGCGTCGTGTA"
LT06_OLIGO_INDEX = "ATTGGC"
RYSHAPE_P7_LT06 = il.P7 + LT06_OLIGO_INDEX + il.TRUSEQ_READ2
P5_GMAC = il.P5[:28] + "TT" + GM_PRIMER + "AC"
READ1_DROPQUARTZ = "ACATT" + GM_PRIMER + "AC"
TPC2 = il.P7[:22]


def rt_primer(barcode_nt: int = 15) -> list[Segment]:
    return [seg("M handle", M, "tso"), seg("AC", "AC"),
            seg("cell barcode", "B" * barcode_nt, "cbc", placeholder=True),
            seg("UMI", "U" * 8, "umi", placeholder=True), seg("dT24", "T" * 24)]


def tagging_primer() -> list[Segment]:
    return [seg("M handle", M, "tso"), seg("dT23", "T" * 23)]


def mrna_rt_scene() -> Scene:
    mrna = [seg("RNA body", "X" * 28, placeholder=True), seg("poly(A)", "A" * 24)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA")
    sc.anneal("RT primer", rt_primer(), to="mRNA", pair=("dT24", "poly(A)"),
              label="barcoded RT primer")
    sc.arrow("RT primer", "reverse transcriptase")
    return sc


def tailed_first_strand() -> Construct:
    return Construct([
        *rt_primer(), seg("cDNA", "X" * 28, placeholder=True),
        seg("TdT poly(A)", "A" * 24),
    ], name="poly(A)-tailed Quartz-Seq2 first strand")


def tagging_scene() -> Scene:
    target = tailed_first_strand()
    sc = Scene()
    sc.strand("first strand", list(target), label="first strand")
    sc.anneal("tagging primer", tagging_primer(), to="first strand",
              pair=("dT23", "TdT poly(A)"), label="tagging primer")
    sc.arrow("tagging primer", "second-strand synthesis")
    return sc


def truncated_adapter_scene() -> Scene:
    """Source-listed truncated Y adapter; Scene enforces its 12-base stem."""
    long = [seg("P7", il.P7, "p7"), seg("i7", LT06_OLIGO_INDEX, "cbc"),
            seg("Read 2 arm", il.TRUSEQ_READ2[:-13], "r2"),
            seg("stem complement", il.TRUSEQ_READ2[-13:-1], "r2"), seg("T", "T")]
    short = [seg("stem", RYSHAPE_P5[:12], "r1"),
             seg("unpaired Read 1 arm", RYSHAPE_P5[12:], "r1")]
    sc = Scene()
    sc.strand("long", long, label="P7 strand")
    sc.anneal("short", short, to="long", pair=("stem", "stem complement"),
              label="short strand", unpaired=("unpaired Read 1 arm",))
    sc.mark("long", "T", "3'-T ligation overhang")
    return sc


READ1_PRIMER = sp.custom("Read 1", "Read1DropQuartz", READ1_DROPQUARTZ,
                         "Sasagawa et al. 2018, Supplementary Table S4")
SEQ_PRIMERS = (READ1_PRIMER, sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])


def final_library(insert_nt: int = 36) -> Construct:
    lib = Construct([
        seg("truncated P5", il.P5[:28], "p5"),
        seg("TT", "TT"),
        seg("gM primer", GM_PRIMER, "tso"),
        seg("AC", "AC"),
        seg("cell barcode", "B" * 15, "cbc", placeholder=True),
        seg("UMI", "U" * 8, "umi", placeholder=True),
        seg("poly(T)", "T" * 24),
        seg("3-prime cDNA fragment", "X" * insert_nt, placeholder=True),
        seg("TruSeq Read 2 reverse complement", revcomp(il.TRUSEQ_READ2), "r2"),
        seg("i7", revcomp(LT06_OLIGO_INDEX), "cbc"),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="Quartz-Seq2 sequencing library")
    for primer in SEQ_PRIMERS:
        if sp.locate(lib, primer) is None:
            raise ValueError(f"{primer.name} does not land on {lib.name}")
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
            sp.feature_spans(final_library(), SEQ_PRIMERS, {"Read 1": 23})]
