"""Molecular model for the MboI implementation of Rao et al. in situ Hi-C."""
from __future__ import annotations
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

import illumina as il
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, complement_segments
from restriction import Digest, MBOI


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    return Segment(name=name, top=top, tag=tag, **kw)


DEMO = Digest(MBOI, "X" * 18, "X" * 18)
FILLED = DEMO.fill_in(biotin_base="A")
SHEAR_RANGE_BP = (300, 500)


def cleavage_scene() -> Scene:
    top = [seg("left locus", DEMO.left_flank, placeholder=True),
           seg("MboI site", MBOI.site, "me"),
           seg("right locus", DEMO.right_flank, placeholder=True)]
    sc = Scene.duplex(top, label="chromatin DNA")
    sc.junction("top", "left locus", "MboI site", "top cut")
    sc.junction("bottom", "MboI site'", "right locus'", "bottom cut")
    sc.labels("top")
    return sc


def cut_left_scene() -> Scene:
    left = seg("left locus", DEMO.left_flank, placeholder=True)
    over = seg("5-prime GATC overhang", MBOI.overhang, "me")
    sc = Scene(); sc.strand("top", [left], label="left fragment")
    sc.anneal("bottom", [over, *complement_segments([left])], to="top",
              pair=("left locus'", "left locus"), label="left fragment")
    sc.mark("bottom", "5-prime GATC overhang", "cohesive end")
    return sc


def cut_right_scene() -> Scene:
    over = seg("5-prime GATC overhang", MBOI.overhang, "me")
    right = seg("right locus", DEMO.right_flank, placeholder=True)
    sc = Scene(); sc.strand("top", [over, right], label="right fragment")
    sc.anneal("bottom", complement_segments([right]), to="top",
              pair=("right locus'", "right locus"), label="right fragment")
    sc.mark("top", "5-prime GATC overhang", "cohesive end")
    return sc


def filled_left_scene() -> Scene:
    made = FILLED.left_new_top
    bio = FILLED.biotin_left_top[0]
    top = [seg("left locus", DEMO.left_flank, placeholder=True),
           seg("new fill before biotin", made[:bio], "me"),
           seg("new biotin-dA", made[bio:bio + 1], "w1"),
           seg("new fill after biotin", made[bio + 1:], "me")]
    sc = Scene.duplex(top, label="filled left end")
    sc.mark("top", "new fill before biotin", "Klenow fill-in",
            through="new fill after biotin")
    sc.mark("top", "new biotin-dA", "biotin")
    return sc


def filled_right_scene() -> Scene:
    # The new strand is stored 5'->3'; reverse its coordinates for the bottom
    # strand as drawn left-to-right (3'->5').
    overhang = MBOI.overhang
    drawn_bio = len(overhang) - 1 - FILLED.biotin_right_bottom[0]
    top = [seg("existing before opposite biotin", overhang[:drawn_bio], "me"),
           seg("opposite new biotin-dA", overhang[drawn_bio:drawn_bio + 1], "w1"),
           seg("existing after opposite biotin", overhang[drawn_bio + 1:], "me"),
           seg("right locus", DEMO.right_flank, placeholder=True)]
    sc = Scene.duplex(top, label="filled right end")
    sc.mark("bottom", "opposite new biotin-dA'", "new bottom-strand biotin-dA")
    return sc


def contact_junction() -> Construct:
    """One chimeric molecule after two filled MboI ends ligate."""
    junction = FILLED.junction
    top_bio = FILLED.junction_biotin_top[0]
    bottom_bio = FILLED.junction_biotin_bottom[0]
    boundary = MBOI.bottom_cut
    return Construct([
        seg("locus A", "X" * 22, placeholder=True),
        seg("left filled prefix", junction[:top_bio], "me"),
        seg("top biotin-dA", junction[top_bio:top_bio + 1], "w1"),
        seg("left filled suffix", junction[top_bio + 1:boundary], "me"),
        seg("right filled prefix", junction[boundary:bottom_bio], "me"),
        seg("opposite bottom biotin-dA", junction[bottom_bio:bottom_bio + 1], "w1"),
        seg("right filled suffix", junction[bottom_bio + 1:], "me"),
        seg("locus B", "X" * 22, placeholder=True),
    ], name="MboI proximity-ligation product")


def sheared_insert() -> Construct:
    return Construct([
        seg("locus A read end", "X" * 28, placeholder=True),
        *list(contact_junction().slice("left filled prefix", "right filled suffix")),
        seg("locus B read end", "X" * 28, placeholder=True),
    ], name="biotin-selected Hi-C insert")


SEQ_PRIMERS = (sp.TRUSEQ["R1"], sp.TRUSEQ["I1"], sp.TRUSEQ["R2"])


def final_library() -> Construct:
    """PCR-completed single-index Illumina library.

    Rao's protocol names an Illumina indexed adapter and Illumina PCR primers without
    printing their sequences. The canonical TruSeq arms are consequently marked inferred.
    """
    lib = Construct([
        seg("P5", il.P5, "p5", inferred=True),
        seg("Read 1 arm", il.TRUSEQ_READ1[4:], "r1", inferred=True),
        *list(sheared_insert()), seg("dA junction", "A", inferred=True),
        seg("Index 1 / Read 2 arm", il.INDEX1_PRIMER, "r2", inferred=True),
        seg("i7 reverse complement", "I" * 8, "cbc", placeholder=True, inferred=True),
        seg("P7 reverse complement", il.P7_RC, "p7", inferred=True),
    ], name="in situ Hi-C sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS,
                         required_roles=("Read 1", "Index 1 (i7)", "Read 2"))
    if problems:
        raise ValueError("invalid in situ Hi-C final library: " + "; ".join(problems))
    return lib
