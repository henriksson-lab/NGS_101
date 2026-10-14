"""Physical duplex-preservation models for rare-variant sequencing."""
from __future__ import annotations

import illumina as il
from chemdraw import Construct, Scene, Segment, feature, oligo, revcomp


# Ding et al. 2023, Supplementary Table 5.  The paper publishes eight matched
# sample-index sets; set 1 is the representative molecule drawn here.  The three Ns
# are the per-molecule UMI, not an undisclosed sequence.
CODEC_SID5_1 = "CTTGAACGGACTGTCCAC"
CODEC_SID7_1 = "CACCGAGCGTTAGACTAC"
CODEC_ADAP5_1 = il.P5 + CODEC_SID5_1 + "T"
CODEC_ADAP7_1 = il.P7 + CODEC_SID7_1 + "T"
CODEC_LINKER_CORE_5 = "GCTAGACTCTGACGTGTTGATCCTCGAAGC"
CODEC_LINKER_CORE_7 = revcomp(CODEC_LINKER_CORE_5)
CODEC_BRID5_1 = (revcomp(CODEC_SID5_1) + "NNN" + revcomp(il.TRUSEQ_READ1)
                 + "TTACATAGTTATCC" + CODEC_LINKER_CORE_5)
CODEC_BRID7_1 = (revcomp(CODEC_SID7_1) + "NNN" + revcomp(il.TRUSEQ_READ2)
                 + "CAATCTATAAGTTGC" + CODEC_LINKER_CORE_7)

# These relationships are the molecular assembly, not regression-test restatements.
# Keeping them beside construction prevents a mistyped published oligo from silently
# becoming a plausible-looking quadruplex.
if not CODEC_BRID5_1.startswith(revcomp(CODEC_SID5_1) + "NNN" + revcomp(il.TRUSEQ_READ1)):
    raise ValueError("CODEC brid5 does not carry the published SID/UMI/Read-1 structure")
if not CODEC_BRID7_1.startswith(revcomp(CODEC_SID7_1) + "NNN" + revcomp(il.TRUSEQ_READ2)):
    raise ValueError("CODEC brid7 does not carry the published SID/UMI/Read-2 structure")
if CODEC_LINKER_CORE_7 != revcomp(CODEC_LINKER_CORE_5):
    raise ValueError("CODEC bridge oligos do not form the published 30-bp linker core")


def codec_oligos():
    """Copyable representative CODEC adapter set, preserving ordering modifications."""
    return (
        oligo("LD4-adap5-1", [Segment("published ordering sequence",
                                      CODEC_ADAP5_1[:-1] + "*T", placeholder=True)],
              three="-3′"),
        oligo("LD4-adap7-1", [Segment("published ordering sequence",
                                      CODEC_ADAP7_1[:-1] + "*T", placeholder=True)],
              three="-3′"),
        oligo("LD4-brid5-1", [Segment("5-prime bridge", CODEC_BRID5_1)],
              five="5′-/5Phos/-", three="-3′"),
        oligo("LD4-brid7-1", [Segment("7-prime bridge", CODEC_BRID7_1)],
              five="5′-/5Phos/-", three="-3′"),
    )


def codec_final_library() -> Construct:
    """Primer-addressable linear view of the linked product in Fig. 1b.

    The source establishes the outward read order on both sides of the central bridge:
    primer, 3-nt UMI, 18-nt matched sample index, terminal T, original-strand copy.
    The bridge itself stays a role token because the figure does not print the extended
    product as one contiguous strand; the four input oligos above remain fully exact.
    """
    sid_set = "CODEC Supplementary Table 5 matched LD4 sample-index sets"
    return Construct([
        Segment("P7", il.P7, "p7"),
        Segment("Crick-strand copy", "x" * 28, placeholder=True),
        Segment("adap7 terminal T", "T"),
        Segment("SID7", CODEC_SID7_1, "cbc",
                feature=feature("sample_index_sid7", "sample_index", "whitelist",
                                whitelist=sid_set)),
        Segment("Read 2 UMI", "N" * 3, "umi",
                feature=feature("umi_r2", "umi", "random")),
        Segment("Read 2 primer site reverse complement", revcomp(il.TRUSEQ_READ2), "r2"),
        Segment("CODEC bridge linker", "[CODEC bridge linker]", "me",
                placeholder=True),
        Segment("Read 1 primer site", il.TRUSEQ_READ1, "r1"),
        Segment("Read 1 UMI", "N" * 3, "umi",
                feature=feature("umi_r1", "umi", "random")),
        Segment("SID5", CODEC_SID5_1, "cbc",
                feature=feature("sample_index_sid5", "sample_index", "whitelist",
                                whitelist=sid_set)),
        Segment("adap5 terminal T", "T"),
        Segment("Watson-strand copy", "X" * 28, placeholder=True),
        Segment("P5 reverse complement", il.P5_RC, "p5"),
    ], name="CODEC linked-duplex sequencing library")


def _p(name: str, base: str, n: int, tag: str | None = None) -> Segment:
    return Segment(name, base * n, tag, placeholder=True)


def codec_linked_duplex() -> Scene:
    """One CODEC product containing readable copies of both original strands."""
    parts = [_p("Watson copy", "X", 28), _p("CODEC linkage", "L", 14, "me"),
             _p("Crick copy", "x", 28)]
    sc = Scene(); sc.strand("linked duplex", parts, label="single CODEC cluster molecule")
    sc.junction("linked duplex", "Watson copy", "CODEC linkage", "strand linkage")
    sc.junction("linked duplex", "CODEC linkage", "Crick copy", "strand linkage")
    sc.note("linked duplex", "one read pair carries both complements of the original duplex")
    sc.labels("linked duplex")
    return sc


def nanoseq_restrictive_end_repair() -> Scene:
    """NanoSeq end state: intact restriction ends, with nick extension suppressed."""
    top = [_p("restriction fragment", "X", 34), Segment("blocked nick", "H", placeholder=True)]
    sc = Scene(); sc.strand("top", top, label="restriction-generated genomic duplex")
    sc.anneal("bottom", [_p("restriction fragment partner", "x", 34)], to="top",
              pair=("restriction fragment partner", "restriction fragment"),
              label="restriction-generated genomic duplex")
    sc.mark("top", "blocked nick", "ddNTP terminates nick extension")
    return sc
