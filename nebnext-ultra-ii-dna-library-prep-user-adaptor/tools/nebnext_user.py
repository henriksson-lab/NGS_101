"""NEBNext Ultra II DNA library prep with the E7600 USER-cleavable adaptor."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Segment, feature, revcomp
from dumbbell import Dumbbell
from endprep import dA_tailed_scene, repair_and_dA_tail


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


INSERT = Construct([seg("DNA insert", "X" * 42, placeholder=True)], name="fragmented DNA")
END_PREP = repair_and_dA_tail(INSERT)
HAIRPIN = il.NEBNEXT_HAIRPIN_MODEL
OPENED = il.NEBNEXT_USER_OPENED
PRE_USER = Dumbbell(INSERT, HAIRPIN.sequence, HAIRPIN.sequence,
                     name="hairpin-ligated NEBNext intermediate")
I5_NAME, I7_NAME = "i501", "i701"
I5 = il.NEBNEXT_I5_SET1[I5_NAME]
I7_IN_OLIGO = il.NEBNEXT_I7_SET1[I7_NAME]
I7_READ = revcomp(I7_IN_OLIGO)
END_PREP_TIMES = ((20, 30), (65, 30))
LIGATION = (20, 15)
USER = (37, 15)
MINIMUM_PCR_CYCLES = 3


def hairpin_oligo() -> Construct:
    p = il.NEBNEXT_HAIRPIN_DU_POS
    n = HAIRPIN.stem_len
    end = len(HAIRPIN.sequence) - HAIRPIN.overhang_3_len
    return Construct([
        seg("5-prime stem", HAIRPIN.sequence[:n], "r2"),
        seg("Read-2-side loop", HAIRPIN.sequence[n:p], "r2"),
        seg("dU", HAIRPIN.sequence[p:p + 1], "w1"),
        seg("Read-1-side loop", HAIRPIN.sequence[p + 1:end - n], "r1"),
        seg("3-prime stem", HAIRPIN.sequence[end - n:end], "r1"),
        seg("3-prime dT", HAIRPIN.overhang_3, None),
    ], name="NEBNext Adaptor for Illumina")


def final_library() -> Construct:
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("i5", I5, "cbc", feature=feature("sample_i5", "sample_index", "whitelist", whitelist="NEBNext index set 1")),
        seg("Read 1 arm", il.TRUSEQ_READ1, "r1"),
        seg("insert", "X" * 42, placeholder=True),
        seg("dA junction", "A"),
        seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2"),
        seg("i7 reverse complement", I7_READ, "cbc", feature=feature("sample_i7", "sample_index", "whitelist", whitelist="NEBNext index set 1")),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="PCR-completed NEBNext library")
    problems = sp.verify(lib, sequencing_primers())
    if problems:
        raise ValueError("invalid NEBNext final library: " + "; ".join(problems))
    return lib


def sequencing_primers() -> tuple[sp.SeqPrimer, ...]:
    return (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["I2"], sp.TRUSEQ["R2"])


def end_prep_scene():
    return dA_tailed_scene(END_PREP)
