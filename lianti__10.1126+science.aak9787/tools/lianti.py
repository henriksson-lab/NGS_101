"""Construct model for Chen et al. LIANTI (Science, 2017).

The defining paper establishes the architecture and reaction order, but the accessible
main text does not print oligo sequences.  LIANTI-specific bases below are transcribed
from the upstream scg_lib_structs page and are therefore kept visibly secondary-source
in the schematic.  Standard Tn5 and Illumina/NEBNext parts come from the shared library.
"""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct, Segment, feature, revcomp

UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "unknown")
FEATURES = {"UMI": UMI_FEATURE, "i7": I7_FEATURE}
from endprep import dA_tailed_scene, repair_and_dA_tail


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


# LIANTI-specific sequence: secondary source (upstream scg_lib_structs), not the paper.
SPACER = "GAACAGAATT"
T7_PROMOTER = "TAATACGACTCACTATA"
T7_START = "GGG"
TRANSPOSON = nx.ME_RC + SPACER + T7_PROMOTER + T7_START + nx.ME
UMI_NT = 8

# Library-prep assignment is also only available in the upstream page.  These bases are
# canonical vendor sequences, imported rather than copied; the upstream specifies 6-nt i7.
I7_NT = 6
NEBNEXT_HAIRPIN = il.NEBNEXT_HAIRPIN_MODEL
NEBNEXT_OPENED = il.NEBNEXT_USER_OPENED


def transposon_segments(*, inferred: bool = True) -> list[Segment]:
    """The single 68-nt hairpin-forming transposon, written 5' to 3'."""
    return [
        seg("ME reverse complement", nx.ME_RC, "me", inferred=inferred),
        seg("spacer", SPACER, inferred=inferred),
        seg("T7 promoter", T7_PROMOTER, "t7", inferred=inferred),
        seg("T7 start", T7_START, "t7", inferred=inferred),
        seg("mosaic end", nx.ME, "me", inferred=inferred),
    ]


def second_strand_primer_segments() -> list[Segment]:
    return [
        seg("UMI", "N" * UMI_NT, "umi", placeholder=True, inferred=True),
        seg("T7 start", T7_START, "t7", inferred=True),
        seg("mosaic end", nx.ME, "me", inferred=True),
    ]


def gap_filled_fragment(insert_nt: int = 28) -> Construct:
    """A representative ds fragment after polymerase gap fill.

    The second end is the reverse complement of the first, so the two T7 promoters point
    inward.  All LIANTI-specific bases stay dotted because their exact sequence is not in
    the accessible defining paper.
    """
    insert = "X" * insert_nt if insert_nt <= 30 else "XXXXXXXXXXXX...XXXXXXXXXXXX"
    right = revcomp(TRANSPOSON)
    lib = Construct([
        *transposon_segments(),
        seg("genomic insert", insert, placeholder=True),
        seg("ME reverse complement, opposite", right[:len(nx.ME)], "me", inferred=True),
        seg("T7 start, opposite", right[len(nx.ME):len(nx.ME) + len(T7_START)],
            "t7", inferred=True),
        seg("T7 promoter, opposite",
            right[len(nx.ME) + len(T7_START):len(nx.ME) + len(T7_START) + len(T7_PROMOTER)],
            "t7", inferred=True),
        seg("spacer, opposite",
            right[len(nx.ME) + len(T7_START) + len(T7_PROMOTER):
                  len(nx.ME) + len(T7_START) + len(T7_PROMOTER) + len(SPACER)],
            inferred=True),
        seg("mosaic end, opposite", right[-len(nx.ME):], "me", inferred=True),
    ], name="gap-filled LIANTI fragment")
    return lib


def transcript(insert_nt: int = 28) -> Construct:
    """One IVT product in DNA letters (T represents U in the actual RNA)."""
    insert = "X" * insert_nt if insert_nt <= 30 else "XXXXXXXXXXXX...XXXXXXXXXXXX"
    return Construct([
        seg("5' GGG", T7_START, "t7", inferred=True),
        seg("5' mosaic end", nx.ME, "me", inferred=True),
        seg("genomic RNA", insert, placeholder=True),
        seg("self-primer site", nx.ME_RC, "me", inferred=True),
        seg("CCC", revcomp(T7_START), "t7", inferred=True),
        seg("promoter reverse complement", revcomp(T7_PROMOTER), "t7", inferred=True),
        seg("spacer reverse complement", revcomp(SPACER), inferred=True),
        seg("3' self-primer", nx.ME, "me", inferred=True),
    ], name="LIANTI RNA transcript")


def umi_amplicon(insert_nt: int = 28) -> Construct:
    insert = "X" * insert_nt if insert_nt <= 30 else "XXXXXXXXXXXX...XXXXXXXXXXXX"
    return Construct([
        *second_strand_primer_segments(),
        seg("genomic copy", insert, placeholder=True),
    ], name="UMI-tagged LIANTI amplicon")


def end_prepared_amplicon(insert_nt: int = 28):
    """The standard NEBNext-polished, dA-tailed LIANTI amplicon."""
    return repair_and_dA_tail(umi_amplicon(insert_nt), inferred=True)


def end_prep_scene(insert_nt: int = 28):
    return dA_tailed_scene(end_prepared_amplicon(insert_nt), label="LIANTI amplicon")


SEQ_PRIMERS = [sp.TRUSEQ[k] for k in ("R1", "I1", "I2", "R2")]


def final_library(insert_nt: int = 28) -> Construct:
    """One of the two possible adapter-ligation orientations, P5 to P7'."""
    insert = "X" * insert_nt if insert_nt <= 30 else "XXXXXXXXXXXX...XXXXXXXXXXXX"
    d_a = end_prepared_amplicon(insert_nt).three_prime_overhang
    lib = Construct([
        # P5 overlaps the first four bases of the Read 1 site (ACAC).
        seg("P5", il.P5[:-4], "p5"),
        seg("Read 1 site", NEBNEXT_OPENED.right_arm, "r1"),
        seg("UMI", "N" * UMI_NT, "umi", placeholder=True, inferred=True),
        seg("T7 start", T7_START, "t7", inferred=True),
        seg("mosaic end", nx.ME, "me", inferred=True),
        seg("genomic insert", insert, placeholder=True),
        seg("dA", d_a, inferred=True),
        seg("index-read site", NEBNEXT_OPENED.left_arm
            + il.INDEX1_PRIMER[len(NEBNEXT_OPENED.left_arm):], "r2"),
        seg("i7", "I" * I7_NT, "cbc", placeholder=True, inferred=True),
        seg("P7'", il.P7_RC, "p7"),
    ], name="LIANTI sequencing library, UMI at P5 end")
    problems = sp.verify(lib, SEQ_PRIMERS)
    if problems:
        raise ValueError("invalid LIANTI final library: " + "; ".join(problems))
    return lib
