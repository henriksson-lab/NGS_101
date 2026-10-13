"""Molecular construct model for Seq-Well (Gierahn et al., 2017)."""
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
I7_FEATURE = feature("sample_index_i7", "sample_index", "unknown")
FEATURES = {"cell barcode": CELL_BARCODE, "UMI": UMI_FEATURE,
            "i7 reverse complement": I7_FEATURE}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


LINKER = "T" * 7
SPACER = "GCCTGTCCGCGG"
BEAD_TEMPLATE = LINKER + rt.SMART_HANDLE + "AC" + "J" * 12 + "N" * 8 + "T" * 30
TSO = rt.SMART_HANDLE + "GAAT" + "GGG"
SMART_PCR = rt.SMART_HANDLE
P5_HYBRID = il.P5 + SPACER + rt.SMART_HANDLE + "AC"
CUSTOM_R1 = SPACER + rt.SMART_HANDLE + "AC"


def bead_primer() -> list[Segment]:
    return [seg("bead linker", LINKER), seg("SMART handle", rt.SMART_HANDLE, "tso"),
            seg("AC", "AC"), seg("cell barcode", "B" * 12, "cbc", placeholder=True),
            seg("UMI", "U" * 8, "umi", placeholder=True), seg("dT30", "T" * 30)]


def tso() -> list[Segment]:
    return [seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("GAAT", "GAAT"),
            seg("rGrGrG", "GGG", "tso")]


def p5_hybrid() -> list[Segment]:
    return [seg("P5", il.P5, "p5"), seg("spacer", SPACER, "r1"),
            seg("SMART handle", rt.SMART_HANDLE, "tso"), seg("AC", "AC")]


def capture_scene() -> Scene:
    mrna=[seg("RNA body","X"*28,placeholder=True),seg("poly(A)","A"*30)]
    sc=Scene(); sc.strand("mRNA",mrna,label="captured mRNA")
    sc.anneal("bead primer",bead_primer(),to="mRNA",pair=("dT30","poly(A)"),label="bead primer",mod5="bead")
    sc.arrow("bead primer","reverse transcriptase"); return sc


def first_strand() -> Construct:
    return Construct([*bead_primer(),seg("antisense cDNA","X"*28,placeholder=True),
                      seg("CCC","CCC","tso")],name="Seq-Well first-strand cDNA")


def template_switch_scene() -> Scene:
    first=first_strand(); sc=Scene(); sc.strand("first strand",list(first),label="first strand",mod5="bead")
    sc.anneal("TSO",tso(),to="first strand",pair=("rGrGrG","CCC"),label="TSO",above=True)
    sc.arrow("first strand","RT switches template"); return sc


def wta_product() -> Construct:
    return Construct([*first_strand(),seg("copied GAAT",revcomp("GAAT")),
                      seg("copied SMART handle",revcomp(rt.SMART_HANDLE),"tso")],
                     name="Seq-Well template-switched cDNA")


def selected_tagmented_strand(insert_nt: int=34) -> Construct:
    return Construct([
        seg("SMART handle",rt.SMART_HANDLE,"tso"),seg("AC","AC"),
        seg("cell barcode","B"*12,"cbc",placeholder=True),seg("UMI","U"*8,"umi",placeholder=True),
        seg("poly(T)","T"*30),seg("antisense cDNA","X"*insert_nt,placeholder=True),
        seg("mosaic end reverse complement",nx.ME_RC,"me"),seg("s7 reverse complement",nx.S7_RC,"s7")],
        name="Seq-Well bead-end tagmented strand")


READ1=sp.custom("Read 1","Seq-Well custom Read 1",CUSTOM_R1,
                "Gierahn et al. 2017, Supplementary Table 1")
SEQ_PRIMERS=(READ1,sp.NEXTERA["I1"],sp.NEXTERA["R2"])


def final_library(insert_nt: int=34) -> Construct:
    selected=selected_tagmented_strand(insert_nt)
    lib=Construct([seg("P5",il.P5,"p5"),seg("spacer",SPACER,"r1"),*list(selected),
                   seg("i7 reverse complement","I"*8,"cbc",placeholder=True),
                   seg("P7 reverse complement",il.P7_RC,"p7")],name="Seq-Well sequencing library")
    for p in SEQ_PRIMERS:
        if sp.locate(lib,p) is None: raise ValueError(f"{p.name} does not land on {lib.name}")
    return lib


def primer_landings(lib: Construct|None=None):
    lib=lib or final_library(); out=[]
    for p in SEQ_PRIMERS:
        h=sp.locate(lib,p)
        if h is None: raise ValueError(f"{p.name} does not land on {lib.name}")
        out.append((p,h))
    return out


def read1_layout() -> list[tuple[str,str]]:
    return [(x.cycles, x.feature.label) for x in
            sp.feature_spans(final_library(), SEQ_PRIMERS, {"Read 1": 20})]


def protocol_graph() -> ProtocolGraph:
    source="Defining source: https://doi.org/10.1038/nmeth.4179"
    graph=ProtocolGraph("Seq-Well",authority=source)
    cells=graph.input("cells","single-cell suspension",MaterialKind.BIOLOGICAL,primary=True)
    beads=graph.input("barcoded_beads","barcoded capture beads",MaterialKind.ADAPTER,
                      state=tuple(bead_primer()))
    membrane=graph.input("membrane","semipermeable microwell membrane",MaterialKind.CONDITION)
    captured=graph.product("captured_rna","microwell-confined bead-captured poly(A) RNA",
                           MaterialKind.MOLECULAR,state=capture_scene())
    empty_wells=graph.product("empty_wells","unproductive cell/bead loading outcomes",
                              MaterialKind.WASTE)
    graph.react("partition_capture","Load cells and beads into membrane-sealed microwells",
                StepKind.PARTITION,
                ((cells,InputRole.TEMPLATE),(beads,InputRole.STRUCTURAL),
                 (membrane,InputRole.CONDITION)),(captured,empty_wells))
    graph.discard(empty_wells,"wells without a productive cell-bead pairing do not yield library")
    reverse_transcriptase=graph.input("reverse_transcriptase","reverse transcriptase",
                                     MaterialKind.CATALYST)
    first=graph.product("first_strand",first_strand().name,MaterialKind.MOLECULAR,
                        state=first_strand())
    graph.react("reverse_transcribe","Synthesize bead-barcoded first-strand cDNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured,InputRole.TEMPLATE),(beads,InputRole.STRUCTURAL),
                 (reverse_transcriptase,InputRole.CATALYST)),
                (first,))
    tso_input=graph.input("tso","template-switching oligo",MaterialKind.OLIGO,
                          state=tuple(tso()))
    switched=graph.product("template_switched",wta_product().name,MaterialKind.MOLECULAR,
                           state=wta_product())
    graph.react("template_switch","Install the distal SMART handle",StepKind.TEMPLATE_SWITCHING,
                ((first,InputRole.TEMPLATE),(tso_input,InputRole.STRUCTURAL),
                 (reverse_transcriptase,InputRole.CATALYST)),(switched,))
    smart=graph.input("smart_primer","SMART PCR primer",MaterialKind.OLIGO,
                      state=(seg("SMART handle",SMART_PCR,"tso"),))
    preamp_pol=graph.input("preamp_polymerase","whole-transcriptome PCR polymerase",
                           MaterialKind.CATALYST)
    amplified=graph.product("amplified_cdna","amplified full-length barcoded cDNA",
                            MaterialKind.MOLECULAR,state=wta_product())
    graph.react("preamp","Amplify template-switched cDNA",StepKind.PCR,
                ((switched,InputRole.TEMPLATE),(smart,InputRole.STRUCTURAL),
                 (preamp_pol,InputRole.CATALYST)),(amplified,),
                geometry=PCRGeometry.symmetric(
                    wta_product(),PCRPrimer("SMART PCR primer",SMART_PCR),
                    statement="Supplementary Table 1 prints the SMART PCR primer and both copied handle sites."))
    transposome=graph.input("transposome","Nextera transposome",MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    tagged=graph.product("selected_tagmented",selected_tagmented_strand().name,
                         MaterialKind.MOLECULAR,state=selected_tagmented_strand())
    tag_waste=graph.product("tagment_waste","tagmentation products without the bead end",
                            MaterialKind.WASTE)
    graph.react("tagment","Tagment cDNA and retain the bead end nearest an s7 insertion",
                StepKind.TAGMENTATION,
                ((amplified,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,tag_waste))
    graph.discard(tag_waste,"fragments lacking the bead-derived handle are not PCR-selected")
    p5=graph.input("p5_primer","P5-SMART selective PCR primer",MaterialKind.OLIGO,
                   state=tuple(p5_hybrid()))
    p7_sequence=il.P7+"N"*8+nx.S7
    p7=graph.input("p7_primer","indexed N7 PCR primer",MaterialKind.OLIGO)
    library_pol=graph.input("library_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final)
    graph.react("library_pcr","Select the bead-end fragment and add P5/P7",
                StepKind.PCR,
                ((tagged,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
                 (p7,InputRole.STRUCTURAL),(library_pol,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.exact(
                    selected_tagmented_strand(),PCRPrimer("P5-SMART",P5_HYBRID,anneal_nt=len(SMART_PCR)+2),
                    PCRPrimer("indexed N7",p7_sequence,anneal_nt=len(nx.S7)),
                    statement="Supplementary Table 1 prints the selective P5-SMART and indexed N7 primer families."))
    lengths={"Read 1":20,"Index 1 (i7)":8,"Read 2":50}
    graph.sequence(library,SequencingPlan(
        library=final,primers=SEQ_PRIMERS,acquired_roles=tuple(lengths),read_lengths=lengths,
        read_uses={"Read 1":"extract","Index 1 (i7)":"demultiplex","Read 2":"include"}))
    graph.declare_reviewed_scope((
        "cells","barcoded_beads","membrane","partition_capture","reverse_transcriptase",
        "reverse_transcribe","tso","template_switch","smart_primer","preamp_polymerase",
        "preamp","transposome","tn5","tagment","p5_primer","p7_primer",
        "library_polymerase","library_pcr","library"),source)
    return graph.require_publishable()
