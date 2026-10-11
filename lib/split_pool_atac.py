"""Exact split-pool ligation chemistry used by SPATAC-seq.

SPATAC-seq installs a two-part first-round barcode with a custom Tn5,
then ligates one further barcode to each end of the same chromatin fragment.
The final molecule is deliberately asymmetric: Read 1 observes round 2 and
the forward half of round 1, whereas Read 2 observes round 3 and the reverse
half of round 1.  This module owns that geometry and validates all four
sequencing-primer sites so pages cannot drift from the demultiplexing layout.

Sequences are from Sun et al. 2024 Supplementary Table 1.  Read coordinates
are independently fixed by the authors' SPATAC-seq processing repository.
"""
from __future__ import annotations

import illumina as il
import nextera as nx
import seqprimers as sp
from batch_ngs import seg
from chemdraw import Construct, MolecularState, Scene, Workflow, feature, revcomp


# Representative oligos, written 5' -> 3'.  Round-1 oligos are phosphorylated;
# the short Tn5 partner is both 5'-phosphorylated and 3'-dideoxy-blocked.
TN5_ME_REV = "CTGTCTCTTATACACAT"
TN5_ME_REV_WRITTEN = "/5Phos/CTGTCTCTTATACACAT/3ddC/"

ROUND1_F_1 = "TATAGCCT"
ROUND1_R_1 = "TTGCTTGC"
ROUND2_1 = "AACGTGAT"
ROUND3_1 = "AACGTGAT"
I5_1 = "ACGGGACT"
I7_1 = "CAGTCCCT"

ROUND1_F_LEADER = "AGGCCAGAGCATTCG"
ROUND1_R_LEADER = "CATCGGCGTACGACT"
ROUND2_PRIMER_TAIL = "ACGACGCTCTTCCGATCT"
ROUND2_BRIDGE = "ATCCACGTGCTTGAG"
ROUND3_PRIMER_SITE = nx.S7 + nx.ME

TN5_F_D1 = ROUND1_F_LEADER + ROUND1_F_1 + nx.ME
TN5_R_D1 = ROUND1_R_LEADER + ROUND1_R_1 + nx.ME
TN5_F_L1 = ROUND2_PRIMER_TAIL + ROUND2_1 + ROUND2_BRIDGE
TN5_F_LA = "TCTGGCCTCTCAAGCACGTGGAT"
TN5_R_L1 = ROUND3_PRIMER_SITE + ROUND3_1
TN5_R_LA_1 = "CGCCGATGATCACGTTCTGTCTC"
I5_PRIMER_1 = il.P5 + I5_1 + il.TRUSEQ_READ1
I7_PRIMER_1 = il.P7 + I7_1 + nx.S7

# The splints prove the two joins rather than merely decorating them.  Their
# 3'->5' display strands must pair across both adjacent oligos.
if revcomp(TN5_F_LA) != ROUND2_BRIDGE + ROUND1_F_LEADER[:8]:
    raise ValueError("SPATAC forward ligation splint does not bridge round 2 to round 1")
if TN5_R_LA_1 != revcomp(ROUND1_R_LEADER[:8]) + revcomp(ROUND3_1) + nx.ME_RC[:7]:
    raise ValueError("SPATAC reverse ligation splint does not bridge round 1 to round 3")


CELL_R1_F = feature(
    "spatac_cell_round1_f", "cell_barcode", "whitelist",
    group="spatac_cell_barcode", part="round 1 forward",
    whitelist="Sun et al. Supplementary Table 1; Lan-lab/SPATAC-seq 48-combination list")
CELL_R1_R = feature(
    "spatac_cell_round1_r", "cell_barcode", "whitelist",
    group="spatac_cell_barcode", part="round 1 reverse",
    whitelist="Sun et al. Supplementary Table 1; Lan-lab/SPATAC-seq 48-combination list")
CELL_R2 = feature(
    "spatac_cell_round2", "cell_barcode", "whitelist",
    group="spatac_cell_barcode", part="round 2",
    whitelist="Sun et al. Supplementary Table 1, Tn5_F_L1-L48")
CELL_R3 = feature(
    "spatac_cell_round3", "cell_barcode", "whitelist",
    group="spatac_cell_barcode", part="round 3",
    whitelist="Sun et al. Supplementary Table 1, Tn5_R_L1-L48")
SAMPLE_I5 = feature(
    "spatac_sample_i5", "sample_index", "whitelist",
    whitelist="Sun et al. Supplementary Table 1, i5 primer panel")
SAMPLE_I7 = feature(
    "spatac_sample_i7", "sample_index", "whitelist",
    whitelist="Sun et al. Supplementary Table 1, i7 primer panel")


def _core(insert_name: str, *, round2: bool, round3: bool) -> list:
    """Top-strand segments from the left barcode arm through the right arm."""
    left = []
    if round2:
        left.extend([
            seg("round 2 primer tail", ROUND2_PRIMER_TAIL, "r1"),
            seg("cell barcode round 2", "B" * 8, "cbc", placeholder=True,
                feature=CELL_R2),
            seg("forward splint-facing bridge", ROUND2_BRIDGE),
        ])
    left.extend([
        seg("round 1 forward leader", ROUND1_F_LEADER),
        seg("cell barcode round 1 forward", "F" * 8, "cbc", placeholder=True,
            feature=CELL_R1_F),
        seg("left mosaic end", nx.ME, "me"),
        seg(insert_name, "X" * 44, placeholder=True),
        seg("right mosaic end reverse complement", nx.ME_RC, "me"),
        seg("cell barcode round 1 reverse complement", "R" * 8, "cbc",
            placeholder=True, feature=CELL_R1_R),
        seg("round 1 reverse leader reverse complement", revcomp(ROUND1_R_LEADER)),
    ])
    if round3:
        left.extend([
            seg("cell barcode round 3 reverse complement", "C" * 8, "cbc",
                placeholder=True, feature=CELL_R3),
            seg("round 3 mosaic end reverse complement", nx.ME_RC, "me"),
            seg("S7 reverse complement", nx.S7_RC, "s7"),
        ])
    return left


def tagmented_fragment(insert_name: str) -> Construct:
    return Construct(_core(insert_name, round2=False, round3=False),
                     name="round-1 barcoded tagmented fragment")


def round2_fragment(insert_name: str) -> Construct:
    return Construct(_core(insert_name, round2=True, round3=False),
                     name="round-2 ligation product")


def round3_fragment(insert_name: str) -> Construct:
    return Construct(_core(insert_name, round2=True, round3=True),
                     name="three-round cell-barcoded fragment")


def library(insert_name: str):
    """PCR-completed hybrid TruSeq/Nextera library with all run sites checked."""
    core = _core(insert_name, round2=True, round3=True)
    # PCR replaces the partial Read-1 priming tail at the extreme left with the
    # complete Read-1 site; the shared 19 bases are not duplicated.
    if core[0].top != ROUND2_PRIMER_TAIL or not il.TRUSEQ_READ1.endswith(core[0].top):
        raise ValueError("SPATAC i5 primer cannot complete the round-2 Read 1 site")
    core = core[1:]
    lib = Construct([
        seg("P5", il.P5, "p5"),
        seg("i5", "J" * 8, "cbc", placeholder=True, feature=SAMPLE_I5),
        seg("TruSeq Read 1 site", il.TRUSEQ_READ1, "r1"),
        *core,
        seg("i7 reverse complement", "I" * 8, "cbc", placeholder=True,
            feature=SAMPLE_I7),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="SPATAC-seq library")
    primers = (
        sp.TRUSEQ["R1"],
        sp.NEXTERA["I1"],
        sp.TRUSEQ["I2"],
        sp.NEXTERA["R2"],
    )
    problems = sp.verify(lib, primers, required_roles=tuple(p.role for p in primers))
    if problems:
        raise ValueError("invalid SPATAC-seq library: " + "; ".join(problems))
    return lib, primers


def _scene(con: Construct, label: str, junctions=()):
    sc = Scene.duplex(list(con), label=label)
    for left, right, text in junctions:
        sc.junction("top", left, right, text)
    sc.labels("top")
    return sc


def tagmented_rows(insert_name: str):
    return _scene(tagmented_fragment(insert_name), "after round-1 indexed Tn5").rows()


def round2_rows(insert_name: str):
    return _scene(round2_fragment(insert_name), "after round-2 ligation", (
        ("forward splint-facing bridge", "round 1 forward leader", "ligation"),
    )).rows()


def round3_rows(insert_name: str):
    return _scene(round3_fragment(insert_name), "after round-3 ligation", (
        ("forward splint-facing bridge", "round 1 forward leader", "ligation"),
        ("round 1 reverse leader reverse complement",
         "cell barcode round 3 reverse complement", "ligation"),
    )).rows()


def final_rows(lib: Construct):
    return _scene(lib, "SPATAC-seq library", (
        ("forward splint-facing bridge", "round 1 forward leader", "ligation"),
        ("round 1 reverse leader reverse complement",
         "cell barcode round 3 reverse complement", "ligation"),
    )).rows()


def workflow(insert_name: str = "accessible genomic DNA") -> Workflow:
    start = Scene.duplex([seg("fixed nuclear chromatin", "X" * 44, placeholder=True)],
                         label="fixed nuclei")
    wf = Workflow(MolecularState("Fixed permeabilized nuclei", tuple(start.rows())))
    wf.react("Split nuclei and transpose with well-indexed custom Tn5",
             tagmented_rows(insert_name), name="Round-1 barcoded chromatin",
             note="Each transferred 5′-phosphorylated adapter carries one 8-nt half of the first cell barcode. The paired halves identify one of 48 used transposomes.")
    wf.react("Pool, split into 48 wells and ligate the round-2 oligo",
             round2_rows(insert_name), name="Round-2 ligation product",
             note="The forward splint anneals across the exact round-2-to-round-1 boundary; ** marks that covalent join.")
    wf.react("Pool, split into 48 wells and ligate the round-3 oligo",
             round3_rows(insert_name), name="Three-round cell-barcoded fragment",
             note="The reverse splint anneals round 3 to the opposite round-1 end. The four 8-nt pieces together form a 32-nt cell identifier.")
    lib, _ = library(insert_name)
    wf.react("Pool, divide into sublibraries, lyse and PCR with indexed primers",
             final_rows(lib), name="Dual-index sequencing library",
             note="PCR completes P5 plus the TruSeq Read-1 site on one end and P7 plus the Nextera S7 site on the other, while adding i5 and i7 sample indexes.")
    return wf


def oligos():
    """Representative public oligos, preserving source modifications."""
    return (
        f"<b>Tn5MErev</b> — 5′-{TN5_ME_REV_WRITTEN}-3′",
        f"<b>Tn5_F_D1</b> — 5′-/5Phos/{TN5_F_D1}-3′",
        f"<b>Tn5_R_D1</b> — 5′-/5Phos/{TN5_R_D1}-3′",
        f"<b>Tn5_F_L1</b> — 5′-{TN5_F_L1}-3′",
        f"<b>Tn5_F_LA</b> — 5′-{TN5_F_LA}-3′",
        f"<b>Tn5_R_L1</b> — 5′-{TN5_R_L1}-3′",
        f"<b>Tn5_R_LA-1</b> — 5′-{TN5_R_LA_1}-3′",
        f"<b>i5 primer-1</b> — 5′-{I5_PRIMER_1}-3′",
        f"<b>i7 primer-1</b> — 5′-{I7_PRIMER_1}-3′",
    )
