"""
Segment and oligo definitions for the Atrandi / PTA single-cell WGS chemistry.

Everything the page draws derives from this file. Sources for each sequence are in the
notes: 01_barcoding_kit.md, 02_library_prep.md, 04_nebnext_illumina.md,
05_barcode_cassette_model.md.

Evidence marking mirrors the notes: segments carrying inferred=True are guesses and render
inside <inf>. Placeholders (barcodes, linkers, insert) are stand-ins, never complemented.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
from chemdraw import Construct, Segment, revcomp

# ------------------------------------------------------- canonical Illumina (verbatim)
# Illumina "Illumina Adapter Sequences", Document # 1000000002694 v22, Sept 2025.

P5_SEQ = "AATGATACGGCGACCACCGAGATCTACAC"
P7_SEQ = "CAAGCAGAAGACGGCATACGAGAT"
TRUSEQ_R1 = "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"
TRUSEQ_R2 = "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"
TRIM_SEEN_IN_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
TRIM_SEEN_IN_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# ------------------------------------------- Atrandi ligation adapter (verbatim, DGPM...206001)
# Top:    /5Phos/GATCGGAAGAGCGTCGTGTAGGGAAAGAGTG*T
# Bottom: /5AmMC6/GCTCTTCCGATCT
LIGADAPT_STEM = "GATCGGAAGAGC"                 # 12 bp duplex stem, top strand
LIGADAPT_ARM = "GTCGTGTAGGGAAAGAGTGT"          # 20 nt single-stranded 3' arm
LIGADAPT_TOP = LIGADAPT_STEM + LIGADAPT_ARM    # 32 nt
LIGADAPT_BOT = "GCTCTTCCGATCT"                 # 13 nt: 12 paired + 1-nt 3'-T overhang

# ------------------------------------------------ Atrandi indexing primers (verbatim)
ATRANDI_I7_INDEX = "AACCTG"                              # 6 nt, i7 only
ATRANDI_I7_ANNEAL = "GTGACTGGAGTTCAGACGTGTGCTCTT"        # 27 nt
ATRANDI_I7 = P7_SEQ + ATRANDI_I7_INDEX + ATRANDI_I7_ANNEAL
ATRANDI_P5_ANNEAL = "ACACTCTTTCCCTACACGAC"              # 20 nt == length of LIGADAPT_ARM
ATRANDI_P5 = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGAC"
ATRANDI_P5_TAIL = ATRANDI_P5[: len(ATRANDI_P5) - len(ATRANDI_P5_ANNEAL)]  # 25 nt 5' overhang

# --------------------------------------------------------------- open design questions
# 05_barcode_cassette_model.md: does the round-D cassette carry the full 34-nt TruSeq
# Read 2 region, or 27 nt (canonical minus the terminal CCGATCT)? Flip this to switch the
# drawing; every diagram re-aligns.
D_ARM_INCLUDES_CCGATCT = False
BARCODE_LEN = 8
LINKER_LEN = 4

_D_ARM_TAIL = TRUSEQ_R2[len(ATRANDI_I7_ANNEAL):]          # "CCGATCT", 7 nt


def _seg(name, top, tag=None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


def read2_arm_segments() -> list[Segment]:
    """The Read-2 / P7 landing site that the round-D cassette must present."""
    segs = [_seg("TruSeq Read 2", ATRANDI_I7_ANNEAL, "t7")]
    if D_ARM_INCLUDES_CCGATCT:
        segs.append(_seg("", _D_ARM_TAIL, "t7", inferred=True,
                         note="present only if the cassette carries the full 34-nt arm"))
    return segs


def barcode_block() -> list[Segment]:
    """D - linker - C - linker - B - linker - A, as read in R2."""
    out: list[Segment] = []
    for i, letter in enumerate("DCBA"):
        out.append(_seg(f"BC-{letter}", letter * BARCODE_LEN, "cbc", placeholder=True))
        if i < 3:
            out.append(_seg("", "L" * LINKER_LEN, "r2", placeholder=True, inferred=True,
                            note="4-nt cohesive overhang, inferred; sequence not published"))
    return out


def final_library() -> Construct:
    """The finished Illumina library, top strand 5'->3' from the P7 end (Atrandi Fig. 2)."""
    segs: list[Segment] = [
        _seg("Illumina P7", P7_SEQ, "p7"),
        _seg("i7", ATRANDI_I7_INDEX, None),
        *read2_arm_segments(),
        *barcode_block(),
        _seg("", "T", None, placeholder=True, inferred=True, bottom="A",
             note="dA/dT ligation junction; Bascet's '+1 to account for ligation'"),
        _seg("insert", "XXX...XXX", None, placeholder=True),
        _seg("TruSeq Read 1", LIGADAPT_TOP, "s5"),
        _seg("Illumina P5", revcomp(ATRANDI_P5_TAIL), "p5"),
    ]
    return Construct(segs, name="Atrandi/PTA scWGS final library")


def r2_read() -> Construct:
    """What Read 2 actually reads: barcode block onward. Offsets must match Bascet."""
    lib = final_library()
    names = [s.name for s in lib]
    start = names.index("BC-D")
    return Construct(lib.segments[start:], name="Read 2")


# ----------------------------------------------------------------- PTA amplicon (06_pta.md)
# Research-grade PTA: exo-resistant random hexamer 5'-NpNpNpNpsNpsN-3' (Dean 2002), extended
# by phi29 and terminated by an alpha-thio-dideoxynucleotide. The commercial kit uses a 6-9mer
# and an undisclosed terminator. 5' end is a hydroxyl, so it needs kinasing before ligation;
# 3' end is a dead end -- no 3'-OH, phosphorothioate linkage, exonuclease-resistant.

PTA_PRIMER_LEN = 6          # 6-9 in the commercial kit (BioSkryb TAS-007 Fig. 2B)
PTA_PRIMER_PS_LINKAGES = 2  # two 3'-terminal phosphorothioates (Thermo SO181 / Dean 2002)


def pta_amplicon() -> Construct:
    """One PTA amplicon as drawn on the page. Both termini are chemically distinctive."""
    return Construct([
        _seg("random primer", "N" * PTA_PRIMER_LEN, "tso", placeholder=True,
             note="5'-OH; last two linkages are phosphorothioate -- needs T4 PNK before ligation"),
        _seg("amplified genomic DNA", "XXX...XXX", None, placeholder=True),
        _seg("terminator", "ddN", "w1", placeholder=True,
             note="alpha-thio-dideoxy: no 3'-OH, PS linkage; cannot be A-tailed or exo-rescued"),
    ], name="PTA amplicon")


def pre_pcr() -> Construct:
    """The barcoded, adapter-ligated molecule as it enters indexing PCR (no P7/i7/P5 yet)."""
    segs = [*read2_arm_segments(), *barcode_block(),
            _seg("", "T", None, placeholder=True, inferred=True, bottom="A"),
            _seg("insert", "XXX...XXX", None, placeholder=True),
            _seg("TruSeq Read 1", LIGADAPT_TOP, "s5")]
    return Construct(segs, name="pre-PCR molecule")


def adapter_both_ends() -> Construct:
    """A fragment that received the ligation adapter at both ends and no barcode."""
    return Construct([
        _seg("Read 1 site", TRUSEQ_R1, "s5"),
        _seg("insert", "XXX...XXX", None, placeholder=True),
        _seg("Read 1 site ", "A" + LIGADAPT_TOP, "s5"),
    ], name="adapter at both ends")
