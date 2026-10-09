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
import seqprimers as sp
from chemdraw import Construct, Segment, feature, revcomp

# ------------------------------------------------------- canonical Illumina (by reference)
# The sequences live once, in lib/illumina.py; these are local aliases, never copies.

P5_SEQ = il.P5
P7_SEQ = il.P7
TRUSEQ_R1 = il.TRUSEQ_READ1
TRUSEQ_R2 = il.TRUSEQ_READ2
TRIM_SEEN_IN_R1 = il.TRIM_SEEN_IN_READ1
TRIM_SEEN_IN_R2 = il.TRIM_SEEN_IN_READ2

# ------------------------------------------- Atrandi ligation adapter (DGPM...206001)
# As printed by Atrandi (citation only -- the values below are derived, and the selftest
# pins the derivation to this text):
#   Top:    /5Phos/GATCGGAAGAGCGTCGTGTAGGGAAAGAGTG*T
#   Bottom: /5AmMC6/GCTCTTCCGATCT
# i.e. the top strand is Illumina's Read-2 trim sequence minus the insert's dA, and the
# bottom strand is the last 13 nt of TruSeq Read 1.
LIGADAPT_TOP = TRIM_SEEN_IN_R2[1:]             # 32 nt
LIGADAPT_STEM = LIGADAPT_TOP[:len(il.STEM)]    # 12 bp duplex stem == il.STEM
LIGADAPT_ARM = LIGADAPT_TOP[len(il.STEM):]     # 20 nt single-stranded 3' arm
LIGADAPT_BOT = TRUSEQ_R1[-13:]                 # 13 nt: 12 paired + 1-nt 3'-T overhang

# ------------------------------------------------ Atrandi indexing primers
# Atrandi's i7 primer = P7 + 6-nt index + the first 27 nt of TruSeq Read 2; its P5 primer =
# P5 + TruSeq Read 1 up to the length of the adapter arm (the two share P5's final ACAC).
ATRANDI_I7_INDEX = "AACCTG"                              # 6 nt, i7 only
ATRANDI_I7_ANNEAL = TRUSEQ_R2[:27]                       # 27 nt
ATRANDI_I7 = P7_SEQ + ATRANDI_I7_INDEX + ATRANDI_I7_ANNEAL
ATRANDI_P5_ANNEAL = TRUSEQ_R1[:len(LIGADAPT_ARM)]        # 20 nt == length of LIGADAPT_ARM
_P5_R1_OVERLAP = 4                                       # "ACAC", end of P5 = start of Read 1
ATRANDI_P5 = P5_SEQ[:-_P5_R1_OVERLAP] + ATRANDI_P5_ANNEAL   # 45 nt
ATRANDI_P5_TAIL = ATRANDI_P5[: len(ATRANDI_P5) - len(ATRANDI_P5_ANNEAL)]  # 25 nt 5' overhang

# ------------------------------------- our indexing primers (as ordered, ref/our_index_primers.tsv)
# Used for both our MDA/PTA scWGS and florian-PTA-rnaseq. Atrandi's i7 design with the 6-nt
# AACCTG swapped for a 10-nt IDT UDP i7 index; the oligo carries the revcomp of the
# sample-sheet i7. P5 is Atrandi's, unchanged. Both end in a 3'-terminal phosphorothioate.
# The selftest checks every modelled oligo against the ordered TSV.
OUR_I7_INDEX = {                                         # index as it sits in the oligo
    "UDP0005": "TAATCTCGTC",
    "UDP0006": "GCGCGATGTT",
    "UDP0007": "AGAGCACTAG",
    "UDP0008": "TGCCTTGATC",
}
OUR_I7_INDEX_LEN = 10
OUR_I7 = {k: P7_SEQ + v + ATRANDI_I7_ANNEAL for k, v in OUR_I7_INDEX.items()}
OUR_P5 = ATRANDI_P5

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
        out.append(_seg(
            f"BC-{letter}", letter * BARCODE_LEN, "cbc", placeholder=True,
            feature=feature(
                f"cell_barcode_{letter.lower()}", "cell_barcode", "combinatorial",
                group="cell_barcode", part=f"round {letter}",
                note="Atrandi split-pool ligation barcode",
            ),
        ))
        if i < 3:
            out.append(_seg("", "L" * LINKER_LEN, "r2", placeholder=True, inferred=True,
                            note="4-nt cohesive overhang, inferred; sequence not published"))
    return out


def dA_junction() -> Segment:
    """The FS-fragment's dA tail opposite the adapter's 3'-T overhang (P5 side). A real base
    pair, drawn explicitly so A + LIGADAPT_TOP == Illumina's Read-2 trim sequence and the
    TruSeq Read 1 primer's 3'-T has its partner (it was formerly hidden in the insert)."""
    return _seg("dA", "A", None, note="insert's dA tail / adapter's 3'-T overhang")


def final_library() -> Construct:
    """The finished Illumina library, top strand 5'->3' from the P7 end (Atrandi Fig. 2)."""
    segs: list[Segment] = [
        _seg("Illumina P7", P7_SEQ, "p7"),
        _seg("i7", ATRANDI_I7_INDEX, None,
             feature=feature("sample_i7", "sample_index", "fixed")),
        *read2_arm_segments(),
        *barcode_block(),
        _seg("", "T", None, placeholder=True, inferred=True, bottom="A",
             note="dA/dT ligation junction; Bascet's '+1 to account for ligation'"),
        _seg("insert", "XXX...XXX", None, placeholder=True),
        dA_junction(),
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
            dA_junction(),
            _seg("TruSeq Read 1", LIGADAPT_TOP, "s5")]
    return Construct(segs, name="pre-PCR molecule")


def adapter_both_ends() -> Construct:
    """A fragment that received the ligation adapter at both ends and no barcode."""
    return Construct([
        _seg("Read 1 site", TRUSEQ_R1, "s5"),
        _seg("insert", "XXX...XXX", None, placeholder=True),
        _seg("Read 1 site ", "A" + LIGADAPT_TOP, "s5"),
    ], name="adapter at both ends")


# ------------------------------------------------------- sequencing primers (lib/seqprimers.py)
def seq_primers() -> list:
    """The primers this library is sequenced with, by reference. Landing sites and first
    bases are computed by seqprimers from final_library(); what depends on the open
    round-D arm question follows D_ARM_INCLUDES_CCGATCT."""
    # Index 1 on the 27-nt arm: its 3' end pairs, its 5' GATCGG is a flap facing barcode D.
    # seqprimers locates and reports that itself, so no declaration is needed either way.
    r2 = sp.TRUSEQ["R2"]
    if not D_ARM_INCLUDES_CCGATCT:
        r2 = sp.mismatching(r2, (
            f"its 3'-terminal {_D_ARM_TAIL} has no partner on the {len(ATRANDI_I7_ANNEAL)}-nt "
            "round-D arm, so it cannot extend. Bascet reads barcode D at Read 2 offset 0, "
            "which implies either the full 34-nt arm (D_ARM_INCLUDES_CCGATCT) or a custom "
            f"Read 2 primer ending ...{ATRANDI_I7_ANNEAL[-6:]}"))
    return [sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["I2"], r2]


def __getattr__(name):
    # SEQ_PRIMERS is recomputed on access so it always follows D_ARM_INCLUDES_CCGATCT.
    if name == "SEQ_PRIMERS":
        return seq_primers()
    raise AttributeError(name)
