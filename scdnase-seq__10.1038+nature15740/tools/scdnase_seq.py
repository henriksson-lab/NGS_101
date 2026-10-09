"""Molecular construct model for scDNase-seq (Pico-Seq).

The paper's cited method identifies Illumina's Multiplexing Sample Prep Oligonucleotide
Kit (ref. 1005709); Illumina publishes that obsolete single-index oligo set.
"""
from __future__ import annotations

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature
from endprep import dA_tailed_scene, repair_and_dA_tail

INDEX_NT = 6


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def digested_fragment_scene() -> Scene:
    """Representative short DNase-I-released genomic fragment."""
    return Scene.duplex([seg("DNase-accessible genomic fragment", "X" * 36,
                             placeholder=True)], label="released DNA")


def a_tailed_fragment_scene() -> Scene:
    """End-repaired fragment with one unpaired 3'-A at each end."""
    core = Construct([seg("end-repaired fragment", "X" * 32, placeholder=True)])
    return dA_tailed_scene(repair_and_dA_tail(core), label="fragment")


def adapter_bottom() -> list[Segment]:
    return [
        seg("Read 1 arm", il.TRUSEQ_READ1[:-(len(il.STEM_COMPLEMENT) + 1)], "r1"),
        seg("stem complement", il.STEM_COMPLEMENT, "r1"),
        seg("3' T", "T", "r1"),
    ]


def adapter_top() -> list[Segment]:
    return [
        seg("stem", il.STEM, "r2"),
        seg("unpaired Read 2 arm", il.INDEX1_PRIMER[len(il.STEM):20], "r2"),
    ]


def adapter_scene() -> Scene:
    """Historical Multiplexing-kit fork; Scene enforces the 12-bp stem."""
    sc = Scene()
    sc.strand("bottom", adapter_bottom(), label="3'-T strand")
    sc.anneal("top", adapter_top(), to="bottom", pair=("stem", "stem complement"),
              label="5'-phosphorylated strand", mod5="p",
              unpaired=("unpaired Read 2 arm",))
    sc.mark("bottom", "3' T", "ligation overhang")
    return sc


def p5_primer() -> list[Segment]:
    return [seg("P5", il.P5, "p5"),
            seg("Read 1 continuation", il.TRUSEQ_READ1[4:], "r1")]


def indexed_p7_primer() -> list[Segment]:
    return [seg("P7", il.P7, "p7"),
            seg("i7", "N" * INDEX_NT, "cbc", placeholder=True),
            seg("Read 2", il.TRUSEQ_READ2, "r2")]


SEQ_PRIMERS = (
    sp.custom("Read 1", "Multiplexing-kit Read 1", il.TRUSEQ_READ1,
              'Illumina "Multiplexing Sample Preparation Oligonucleotide Kit" ref. 1005709'),
    sp.custom("Index 1 (i7)", "Multiplexing-kit index primer", il.INDEX1_PRIMER,
              'Illumina "Multiplexing Sample Preparation Oligonucleotide Kit" ref. 1005709'),
    sp.custom("Read 2", "Multiplexing-kit Read 2", il.TRUSEQ_READ2,
              'Illumina "Multiplexing Sample Preparation Oligonucleotide Kit" ref. 1005709'),
)


def final_library(insert_nt: int = 36) -> Construct:
    """PCR-completed single-index library; required primer sites are validated here."""
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("Read 1 arm", il.TRUSEQ_READ1[4:], "r1"),
        seg("DNase-accessible genomic insert", "X" * insert_nt, placeholder=True),
        seg("dA junction", "A"),
        seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        seg("i7 reverse complement", "N" * INDEX_NT, "cbc", placeholder=True,
            feature=feature("cell_i7", "cell_barcode", "whitelist", whitelist="published indexed primer set")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="scDNase-seq library")
    problems = sp.verify(lib, SEQ_PRIMERS,
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if problems:
        raise ValueError("invalid scDNase-seq library: " + "; ".join(problems))
    return lib


def primer_landings():
    lib = final_library()
    return [(primer, sp.locate(lib, primer)) for primer in SEQ_PRIMERS]
