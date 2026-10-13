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
from adapter_ligation import reported_adapter_end
from pcr import PCRGeometry
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)


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


def protocol_graph(insert_name: str, *, authority: str) -> ProtocolGraph:
    graph=ProtocolGraph("SPATAC-seq",authority=authority)
    nuclei=graph.input("nuclei","fixed permeabilized nuclei",MaterialKind.BIOLOGICAL,primary=True)
    r1a=graph.product("round1_well_a","representative round-1 well A nuclei",MaterialKind.BIOLOGICAL)
    r1b=graph.product("round1_well_b","representative round-1 well B nuclei",MaterialKind.BIOLOGICAL)
    graph.react("round1_split","Split nuclei across custom indexed-Tn5 wells",StepKind.SPLIT,
                ((nuclei,InputRole.TEMPLATE),),(r1a,r1b))
    transposomes=graph.input("round1_transposomes","two-ended round-1 barcode transposome set",
                             MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged_a=graph.product("tagged_a","round-1-barcoded nuclei from well A",MaterialKind.MOLECULAR,
                           state=tagmented_fragment(insert_name))
    tagged_b=graph.product("tagged_b","round-1-barcoded nuclei from well B",MaterialKind.MOLECULAR,
                           state=tagmented_fragment(insert_name))
    waste_a=graph.product("tagment_waste_a","non-library products from well A",MaterialKind.WASTE)
    waste_b=graph.product("tagment_waste_b","non-library products from well B",MaterialKind.WASTE)
    graph.react("tagment_a","Transpose with representative indexed Tn5 well A",StepKind.TAGMENTATION,
                ((r1a,InputRole.TEMPLATE),(transposomes,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),
                (tagged_a,waste_a))
    graph.react("tagment_b","Transpose with representative indexed Tn5 well B",StepKind.TAGMENTATION,
                ((r1b,InputRole.TEMPLATE),(transposomes,InputRole.STRUCTURAL),(tn5,InputRole.CATALYST)),
                (tagged_b,waste_b))
    graph.discard(waste_a,"unusable transposition products do not reach library PCR")
    graph.discard(waste_b,"unusable transposition products do not reach library PCR")
    pool1=graph.product("round1_pool","pooled round-1-barcoded nuclei",MaterialKind.MOLECULAR)
    graph.react("round1_pooling","Pool round-1 tagmentation wells",StepKind.POOL,
                ((tagged_a,InputRole.TEMPLATE),(tagged_b,InputRole.TEMPLATE)),(pool1,))
    r2a=graph.product("round2_well_a","representative round-2 well A nuclei",MaterialKind.MOLECULAR)
    r2b=graph.product("round2_well_b","representative round-2 well B nuclei",MaterialKind.MOLECULAR)
    graph.react("round2_split","Redistribute nuclei across round-2 ligation wells",StepKind.SPLIT,
                ((pool1,InputRole.TEMPLATE),),(r2a,r2b))
    oligo2=graph.input("round2_oligo","round-2 barcode oligo and splint",MaterialKind.ADAPTER)
    ligase=graph.input("ligase","DNA ligase",MaterialKind.CATALYST)
    geometry2=reported_adapter_end(
        "round-2 barcode oligo",kind="splinted",polarity="5-prime",
        terminal_sequence=ROUND2_BRIDGE,
        statement="The published forward splint reverse-complements the round-2 bridge plus the round-1 leader across the ligation bond.")
    out2a=graph.product("round2_product_a","round-2-ligated nuclei from well A",MaterialKind.MOLECULAR,
                        state=round2_fragment(insert_name))
    out2b=graph.product("round2_product_b","round-2-ligated nuclei from well B",MaterialKind.MOLECULAR,
                        state=round2_fragment(insert_name))
    for suffix,template,out in (("a",r2a,out2a),("b",r2b,out2b)):
        graph.react("round2_ligation_"+suffix,"Splint-ligate round-2 barcode in well "+suffix.upper(),
                    StepKind.LIGATION,
                    ((template,InputRole.TEMPLATE),(oligo2,InputRole.STRUCTURAL),
                     (ligase,InputRole.CATALYST)),(out,),geometry=geometry2)
    pool2=graph.product("round2_pool","pooled two-round-barcoded nuclei",MaterialKind.MOLECULAR)
    graph.react("round2_pooling","Pool round-2 ligation wells",StepKind.POOL,
                ((out2a,InputRole.TEMPLATE),(out2b,InputRole.TEMPLATE)),(pool2,))
    r3a=graph.product("round3_well_a","representative round-3 well A nuclei",MaterialKind.MOLECULAR)
    r3b=graph.product("round3_well_b","representative round-3 well B nuclei",MaterialKind.MOLECULAR)
    graph.react("round3_split","Redistribute nuclei across round-3 ligation wells",StepKind.SPLIT,
                ((pool2,InputRole.TEMPLATE),),(r3a,r3b))
    oligo3=graph.input("round3_oligo","round-3 barcode oligo and splint",MaterialKind.ADAPTER)
    geometry3=reported_adapter_end(
        "round-3 barcode oligo",kind="splinted",polarity="5-prime",terminal_sequence=ROUND3_1,
        statement="The published reverse splint spans the round-1 reverse leader, round-3 barcode and mosaic-end sequence.")
    out3a=graph.product("round3_product_a","three-round-barcoded nuclei from well A",MaterialKind.MOLECULAR,
                        state=round3_fragment(insert_name))
    out3b=graph.product("round3_product_b","three-round-barcoded nuclei from well B",MaterialKind.MOLECULAR,
                        state=round3_fragment(insert_name))
    for suffix,template,out in (("a",r3a,out3a),("b",r3b,out3b)):
        graph.react("round3_ligation_"+suffix,"Splint-ligate round-3 barcode in well "+suffix.upper(),
                    StepKind.LIGATION,
                    ((template,InputRole.TEMPLATE),(oligo3,InputRole.STRUCTURAL),
                     (ligase,InputRole.CATALYST)),(out,),geometry=geometry3)
    pool3=graph.product("round3_pool","pooled three-round-barcoded nuclei",MaterialKind.MOLECULAR)
    graph.react("round3_pooling","Pool round-3 ligation wells",StepKind.POOL,
                ((out3a,InputRole.TEMPLATE),(out3b,InputRole.TEMPLATE)),(pool3,))
    sub_a=graph.product("sublibrary_a","PCR sublibrary A input",MaterialKind.MOLECULAR)
    sub_b=graph.product("sublibrary_b","PCR sublibrary B input",MaterialKind.MOLECULAR)
    graph.react("sublibrary_split","Divide pooled nuclei into indexed PCR sublibraries",StepKind.SPLIT,
                ((pool3,InputRole.TEMPLATE),),(sub_a,sub_b))
    i5=graph.input("i5_primer","indexed i5 PCR primer",MaterialKind.OLIGO)
    i7=graph.input("i7_primer","indexed i7 PCR primer",MaterialKind.OLIGO)
    polymerase=graph.input("pcr_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    final,primers=library(insert_name)
    libraries=[]
    for suffix,template in (("a",sub_a),("b",sub_b)):
        out=graph.product("library_"+suffix,final.name+" sublibrary "+suffix.upper(),
                          MaterialKind.LIBRARY,state=final)
        graph.react("library_pcr_"+suffix,"Lyse and PCR-complete sublibrary "+suffix.upper(),
                    StepKind.PCR,
                    ((template,InputRole.TEMPLATE),(i5,InputRole.STRUCTURAL),
                     (i7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(out,),
                    geometry=PCRGeometry.unavailable(final,
                        reason="The exact three-round intermediate is known; PCR completes the verified outer arms."))
        graph.sequence(out,SequencingPlan(
            library=final,primers=primers,acquired_roles=tuple(p.role for p in primers),
            read_lengths={"Read 1":115,"Index 1 (i7)":8,"Index 2 (i5)":8,"Read 2":100},
            read_uses={"Read 1":"extract_then_include","Index 1 (i7)":"demultiplex",
                       "Index 2 (i5)":"demultiplex","Read 2":"extract_then_include"}))
        libraries.append(out)
    graph.declare_reviewed_scope((
        "nuclei","round1_split","round1_transposomes","tn5","tagment_a","tagment_b",
        "round1_pooling","round2_split","round2_oligo","ligase","round2_ligation_a",
        "round2_ligation_b","round2_pooling","round3_split","round3_oligo",
        "round3_ligation_a","round3_ligation_b","round3_pooling","sublibrary_split",
        "i5_primer","i7_primer","pcr_polymerase","library_pcr_a","library_pcr_b",
        "library_a","library_b"),authority)
    return graph.require_publishable()
