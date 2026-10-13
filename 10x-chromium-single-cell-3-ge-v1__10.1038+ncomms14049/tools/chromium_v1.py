"""Architecture model for 10x Chromium Single Cell 3' Gene Expression v1."""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from tenx_graph import chromium_rna_graph

CELL_BARCODE = feature("cell_barcode", "cell_barcode", "whitelist",
                       whitelist="10x-chromium-3prime-v1")
UMI = feature("umi", "umi", "random")
I5 = feature("sample_index_i5", "sample_index", "unknown")
FEATURES = {"cell barcode": CELL_BARCODE, "cell barcode reverse complement": CELL_BARCODE,
            "UMI": UMI, "i5 sample index": I5}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


def inferred(name: str, sequence: str, tag: str | None = None, **kw) -> Segment:
    """The v1 guide publishes architecture but not oligo sequences."""
    return seg(name, sequence, tag, inferred=True, **kw)


def bead_oligo() -> Construct:
    return Construct([
        inferred("P7", il.P7, "p7"),
        seg("cell barcode", "B" * 14, "cbc", placeholder=True),
        inferred("TruSeq Read 2", il.TRUSEQ_READ2, "r2"),
        seg("UMI", "N" * 10, "umi", placeholder=True), seg("poly(T)", "T" * 30),
    ], name="10x 3-prime v1 bead oligo architecture")


def rt_scene(insert_nt: int = 28) -> Scene:
    mrna = [seg("transcript", "X" * insert_nt, placeholder=True), seg("poly(A)", "A" * 30)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA", mod3="poly(A)")
    sc.anneal("bead", list(bead_oligo()), to="mRNA", pair=("poly(T)", "poly(A)"),
              label="released gel-bead primer")
    sc.arrow("bead", "reverse transcription and template switching")
    return sc


def first_strand(insert_nt: int = 28) -> Construct:
    return Construct([*bead_oligo(), seg("antisense cDNA", "X" * insert_nt,
                                        placeholder=True),
                      inferred("template-switch end", "N" * 18, "tso",
                               placeholder=True)],
                     name="10x 3-prime v1 first-strand architecture")


def amplified_cdna(insert_nt: int = 28) -> Construct:
    return Construct(list(first_strand(insert_nt)),
                     name="10x 3-prime v1 amplified-cDNA architecture")


SEQ_PRIMERS = tuple(sp.TRUSEQ[k] for k in ("R1", "I1", "I2", "R2"))


def final_library(insert_nt: int = 28) -> Construct:
    lib = Construct([
        inferred("P5", il.P5, "p5"), seg("i5 sample index", "I" * 8, "cbc", placeholder=True),
        inferred("TruSeq Read 1", il.TRUSEQ_READ1, "r1"),
        seg("3-prime cDNA", "X" * insert_nt, placeholder=True),
        seg("poly(A)", "A" * 30), seg("UMI", "N" * 10, "umi", placeholder=True),
        inferred("TruSeq Read-2 reverse complement", revcomp(il.TRUSEQ_READ2), "r2"),
        seg("cell barcode reverse complement", "B" * 14, "cbc", placeholder=True),
        inferred("P7 reverse complement", il.P7_RC, "p7"),
    ], name="10x Chromium 3-prime v1 library")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid 10x 3-prime v1 library: " + "; ".join(problems))
    return lib


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")


def protocol_graph():
    authority = ('Defining source: <a href="https://doi.org/10.1038/ncomms14049">'
                 'Zheng et al., <i>Nature Communications</i> (2017)</a>; protocol '
                 'authority: 10x Genomics CG00026 Rev B.')
    return chromium_rna_graph(
        title="10x Chromium Single Cell 3-prime Gene Expression v1",
        authority=authority, final=final_library(), primers=SEQ_PRIMERS,
        first_strand=first_strand(), amplified_cdna=amplified_cdna(),
        read_lengths={"Read 1": 98, "Index 1 (i7)": 14,
                      "Index 2 (i5)": 8, "Read 2": 10},
        inferred_intermediates=True)
