"""Molecular construct model for Seq-Well S3 (Hughes et al., 2020)."""
from __future__ import annotations

import illumina as il
import nextera as nx
import rt
import seqprimers as sp
from chemdraw import Construct, Scene, Segment, feature, revcomp
from pcr import PCRGeometry, PCRPrimer
from protocol_graph import (InputRole, MaterialKind, ProtocolGraph, SequencingPlan,
                            StepKind)

CELL_BARCODE = feature("cell_barcode", "cell_barcode", "random")
UMI_FEATURE = feature("umi", "umi", "random")
I7_FEATURE = feature("sample_index_i7", "sample_index", "fixed")
FEATURES = {"cell barcode": CELL_BARCODE, "UMI": UMI_FEATURE,
            "i7 index read": I7_FEATURE}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


BEAD_PREFIX = "T" * 7
BEAD_SUFFIX = "AC"
P5_SPACER = "GCCTGTCCGCGG"
S3_RANDOM_TAIL = "NNNGGNNNB"
I7_INDEX_READ = "TAAGGCGA"
TSO = rt.SMART_HANDLE + "GAAT" + "GGG"
P5_HYBRID = il.P5 + P5_SPACER + rt.SMART_HANDLE + BEAD_SUFFIX
CUSTOM_R1 = P5_SPACER + rt.SMART_HANDLE + BEAD_SUFFIX


def bead_oligo() -> Construct:
    return Construct([
        seg("bead prefix", BEAD_PREFIX), seg("SMART handle", rt.SMART_HANDLE, "tso"),
        seg("bead constant", BEAD_SUFFIX),
        seg("cell barcode", "B" * 12, "cbc", placeholder=True),
        seg("UMI", "N" * 8, "umi", placeholder=True), seg("poly(T)", "T" * 30),
    ], name="Seq-Well S3 bead oligo")


def tso_oligo() -> list[Segment]:
    return [seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("GAAT", "GAAT"),
            seg("rGrGrG", "rGrGrG", placeholder=True)]


def randomer() -> list[Segment]:
    return [seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("GA", "GA"),
            seg("random priming tail", S3_RANDOM_TAIL, placeholder=True)]


def rt_scene(insert_nt: int = 28) -> Scene:
    mrna = [seg("transcript", "X" * insert_nt, placeholder=True),
            seg("poly(A)", "A" * 30)]
    sc = Scene()
    sc.strand("mRNA", mrna, label="mRNA", mod3="poly(A)")
    sc.anneal("bead", list(bead_oligo()), to="mRNA", pair=("poly(T)", "poly(A)"),
              label="bead primer")
    sc.arrow("bead", "reverse transcription; template switching may occur")
    return sc


def random_second_strand_scene() -> Scene:
    # Placeholder text is deliberately symmetric under the Scene's reverse placement:
    # the exact random landing sequence is molecule-specific and is not asserted here.
    template = [seg("first-strand cDNA", "X" * 20, placeholder=True),
                seg("randomer landing", "BNNNCCNNN", placeholder=True)]
    sc = Scene()
    sc.strand("first", template, label="bead-bound first strand")
    sc.anneal("randomer", randomer(), to="first",
              pair=("random priming tail", "randomer landing"), label="S3 randomer")
    sc.arrow("randomer", "Klenow exo-minus makes the second strand toward the bead")
    return sc


def wta_product(insert_nt: int = 28) -> Construct:
    return Construct([
        seg("SMART handle",rt.SMART_HANDLE,"tso"),seg("bead constant",BEAD_SUFFIX),
        seg("cell barcode","B"*12,"cbc",placeholder=True),
        seg("UMI","N"*8,"umi",placeholder=True),seg("poly(T)","T"*30),
        seg("cDNA","X"*insert_nt,placeholder=True),
        seg("copied GA",revcomp("GA")),
        seg("copied SMART handle",revcomp(rt.SMART_HANDLE),"tso")],
        name="Seq-Well S3 handle-flanked WTA product")


def selected_tagmented_strand(insert_nt: int = 28) -> Construct:
    return Construct([
        seg("SMART handle",rt.SMART_HANDLE,"tso"),seg("bead constant",BEAD_SUFFIX),
        seg("cell barcode","B"*12,"cbc",placeholder=True),
        seg("UMI","N"*8,"umi",placeholder=True),seg("poly(T)","T"*30),
        seg("cDNA","X"*insert_nt,placeholder=True),
        seg("mosaic end reverse complement",nx.ME_RC,"me"),
        seg("s7 reverse complement",nx.S7_RC,"s7")],
        name="Seq-Well S3 selected bead-end fragment")


SEQ_PRIMERS = (
    sp.custom("Read 1", "Seq-Well S3 custom Read 1", CUSTOM_R1,
              "Hughes et al. 2020 Key Resources Table"),
    sp.NEXTERA["I1"], sp.NEXTERA["R2"],
)
RUN_ROLES = ("Read 1", "Index 1 (i7)", "Read 2")


def final_library(insert_nt: int = 28) -> Construct:
    lib = Construct([
        seg("P5", il.P5, "p5"), seg("custom spacer", P5_SPACER),
        seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("bead constant", BEAD_SUFFIX),
        seg("cell barcode", "B" * 12, "cbc", placeholder=True),
        seg("UMI", "N" * 8, "umi", placeholder=True), seg("poly(T)", "T" * 30),
        seg("cDNA", "X" * insert_nt, placeholder=True),
        seg("mosaic end reverse complement", nx.ME_RC, "me"),
        seg("s7 reverse complement", nx.S7_RC, "s7"),
        seg("i7 index read", I7_INDEX_READ, "cbc"),
        seg("P7 reverse complement", il.P7_RC, "p7"),
    ], name="Seq-Well S3 sequencing library")
    problems = sp.verify(lib, SEQ_PRIMERS, required_roles=RUN_ROLES)
    if problems:
        raise ValueError("invalid final Seq-Well S3 construct: " + "; ".join(problems))
    return lib


def final_scene() -> Scene:
    return Scene.duplex(list(final_library()), label="library")


def protocol_graph() -> ProtocolGraph:
    source="Defining source: https://doi.org/10.1016/j.immuni.2020.09.015"
    graph=ProtocolGraph("Seq-Well S3",authority=source)
    cells=graph.input("cells","single-cell suspension",MaterialKind.BIOLOGICAL,primary=True)
    beads=graph.input("barcoded_beads","barcoded capture beads",MaterialKind.ADAPTER,
                      state=bead_oligo())
    membrane=graph.input("membrane","semipermeable nanowell membrane",MaterialKind.CONDITION)
    captured=graph.product("captured_rna","nanowell-confined bead-captured poly(A) RNA",
                           MaterialKind.MOLECULAR,state=rt_scene())
    empty=graph.product("empty_wells","unproductive cell/bead loading outcomes",
                        MaterialKind.WASTE)
    graph.react("partition_capture","Load cells and beads into membrane-sealed nanowells",
                StepKind.PARTITION,
                ((cells,InputRole.TEMPLATE),(beads,InputRole.STRUCTURAL),
                 (membrane,InputRole.CONDITION)),(captured,empty))
    graph.discard(empty,"wells without a productive cell-bead pairing do not yield library")
    rt_enzyme=graph.input("reverse_transcriptase","reverse transcriptase",MaterialKind.CATALYST)
    first=graph.product("first_strand","bead-barcoded first-strand cDNA",
                        MaterialKind.MOLECULAR,state=bead_oligo())
    graph.react("reverse_transcribe","Synthesize bead-barcoded first-strand cDNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured,InputRole.TEMPLATE),(beads,InputRole.STRUCTURAL),
                 (rt_enzyme,InputRole.CATALYST)),(first,))
    tso=graph.input("tso","template-switching oligo",MaterialKind.OLIGO,
                    state=tuple(tso_oligo()))
    switched=graph.product("switched_cdna","template-switched cDNA molecules",
                           MaterialKind.MOLECULAR,state=wta_product())
    graph.react("template_switch","Generate the template-switched cDNA fraction",
                StepKind.TEMPLATE_SWITCHING,
                ((first,InputRole.TEMPLATE),(tso,InputRole.STRUCTURAL),
                 (rt_enzyme,InputRole.CATALYST)),(switched,))
    alkali=graph.input("alkali","RNA-stripping alkali",MaterialKind.CONDITION)
    stripped=graph.product("stripped_first_strand","RNA-stripped bead-bound first strand",
                           MaterialKind.MOLECULAR)
    removed_rna=graph.product("removed_rna","hydrolysed RNA template",MaterialKind.WASTE)
    graph.react("strip_rna","Remove RNA from the first-strand cDNA",StepKind.STRAND_SEPARATION,
                ((first,InputRole.TEMPLATE),(alkali,InputRole.CONDITION)),
                (stripped,removed_rna))
    graph.discard(removed_rna,"the RNA template is removed before random second-strand synthesis")
    randomer_input=graph.input("s3_randomer","S3 random second-strand primer",
                               MaterialKind.OLIGO,state=tuple(randomer()))
    klenow=graph.input("klenow","Klenow exo-minus polymerase",MaterialKind.CATALYST)
    rescued=graph.product("random_primed_cdna","random-primer-rescued cDNA molecules",
                          MaterialKind.MOLECULAR,state=wta_product())
    graph.react("random_second_strand","Make a handle-tagged second strand from unswitched cDNA",
                StepKind.EXTENSION,
                ((stripped,InputRole.TEMPLATE),(randomer_input,InputRole.STRUCTURAL),
                 (klenow,InputRole.CATALYST)),(rescued,))
    pooled=graph.product("wta_pool",wta_product().name,MaterialKind.MOLECULAR,
                         state=wta_product())
    graph.react("pool_routes","Pool template-switched and random-primer-rescued products",
                StepKind.POOL,
                ((switched,InputRole.TEMPLATE),(rescued,InputRole.TEMPLATE)),(pooled,),
                changes_molecule=False)
    smart=graph.input("smart_primer","SMART PCR primer",MaterialKind.OLIGO)
    preamp_pol=graph.input("preamp_polymerase","WTA PCR polymerase",MaterialKind.CATALYST)
    amplified=graph.product("amplified_cdna","amplified S3 WTA cDNA",MaterialKind.MOLECULAR,
                            state=wta_product())
    graph.react("preamp","Amplify handle-flanked cDNA",StepKind.PCR,
                ((pooled,InputRole.TEMPLATE),(smart,InputRole.STRUCTURAL),
                 (preamp_pol,InputRole.CATALYST)),(amplified,),
                geometry=PCRGeometry.symmetric(
                    wta_product(),PCRPrimer("SMART PCR primer",rt.SMART_HANDLE),
                    statement="The Key Resources Table prints the SMART primer used at both handle boundaries."))
    transposome=graph.input("transposome","Nextera transposome",MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    selected=selected_tagmented_strand()
    tagged=graph.product("selected_tagmented",selected.name,MaterialKind.MOLECULAR,state=selected)
    waste=graph.product("tagment_waste","tagmentation products without the bead end",
                        MaterialKind.WASTE)
    graph.react("tagment","Tagment amplified cDNA",StepKind.TAGMENTATION,
                ((amplified,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,waste))
    graph.discard(waste,"fragments lacking the bead-derived handle are not PCR-selected")
    p5=graph.input("p5_primer","P5-SMART selective primer",MaterialKind.OLIGO)
    p7=graph.input("p7_primer","N7xx indexed primer",MaterialKind.OLIGO)
    library_pol=graph.input("library_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final)
    graph.react("library_pcr","Select the bead-end/s7 fragment and add P5/P7",StepKind.PCR,
                ((tagged,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
                 (p7,InputRole.STRUCTURAL),(library_pol,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.exact(
                    selected,PCRPrimer("P5-SMART",P5_HYBRID,
                                       anneal_nt=len(rt.SMART_HANDLE)+len(BEAD_SUFFIX)),
                    PCRPrimer("N7xx",il.P7+I7_INDEX_READ+nx.S7,anneal_nt=len(nx.S7)),
                    statement="The Key Resources Table defines the P5-SMART and N7xx primer boundaries."))
    lengths={"Read 1":20,"Index 1 (i7)":8,"Read 2":50}
    graph.sequence(library,SequencingPlan(
        library=final,primers=SEQ_PRIMERS,acquired_roles=RUN_ROLES,read_lengths=lengths,
        read_uses={"Read 1":"extract","Index 1 (i7)":"demultiplex","Read 2":"include"}))
    graph.declare_reviewed_scope((
        "cells","barcoded_beads","membrane","partition_capture","reverse_transcriptase",
        "reverse_transcribe","tso","template_switch","alkali","strip_rna","s3_randomer",
        "klenow","random_second_strand","pool_routes","smart_primer","preamp_polymerase",
        "preamp","transposome","tn5","tagment","p5_primer","p7_primer",
        "library_polymerase","library_pcr","library"),source)
    return graph.require_publishable()
