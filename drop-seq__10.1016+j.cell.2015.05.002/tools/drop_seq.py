"""Molecular model for Drop-seq, with unavailable oligo-table regions inferred by default."""
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
READ_LENGTHS = {"Read 1": 20, "Index 1 (i7)": 8, "Read 2": 50}


def seg(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, **kw)


def inferred(name: str, top: str, tag: str | None = None, **kw) -> Segment:
    """Unavailable Drop-seq Table S6 sequence: always rendered as inferred."""
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, inferred=True, **kw)


LINKER="T"*7
SPACER="GCCTGTCCGCGG"


def bead_primer(batch: str="B") -> list[Segment]:
    if batch not in ("A","B"): raise ValueError("Drop-seq bead batch must be A or B")
    constant="ACGT" if batch=="A" else "AC"
    return [inferred("bead linker",LINKER),inferred("SMART handle",rt.SMART_HANDLE,"tso"),
            inferred(f"batch {batch} constant",constant),
            inferred("cell barcode","B"*12,"cbc",placeholder=True),
            inferred("UMI","U"*8,"umi",placeholder=True),inferred("dT30","T"*30)]


def tso() -> list[Segment]:
    return [inferred("SMART handle",rt.SMART_HANDLE,"tso"),inferred("GAAT","GAAT"),
            inferred("rGrGrG","GGG","tso")]


def p5_hybrid(batch: str="B") -> list[Segment]:
    constant="ACGT" if batch=="A" else "AC"
    return [seg("P5",il.P5,"p5"),inferred("spacer",SPACER,"r1"),
            inferred("SMART handle",rt.SMART_HANDLE,"tso"),inferred(f"batch {batch} constant",constant)]


def custom_r1(batch: str="B") -> str:
    constant="ACGT" if batch=="A" else "AC"
    return SPACER+rt.SMART_HANDLE+constant


def capture_scene(batch: str="B") -> Scene:
    mrna=[seg("RNA body","X"*28,placeholder=True),seg("poly(A)","A"*30)]
    sc=Scene(); sc.strand("mRNA",mrna,label="captured mRNA")
    sc.anneal("bead primer",bead_primer(batch),to="mRNA",pair=("dT30","poly(A)"),label=f"batch {batch} bead",mod5="bead")
    sc.arrow("bead primer","reverse transcriptase"); return sc


def first_strand(batch: str="B") -> Construct:
    return Construct([*bead_primer(batch),seg("antisense cDNA","X"*28,placeholder=True),
                      inferred("CCC","CCC","tso")],name=f"Drop-seq batch {batch} first strand")


def template_switch_scene(batch: str="B") -> Scene:
    first=first_strand(batch); sc=Scene(); sc.strand("first strand",list(first),label="first strand",mod5="bead")
    sc.anneal("TSO",tso(),to="first strand",pair=("rGrGrG","CCC"),label="TSO",above=True)
    sc.arrow("first strand","RT switches template"); return sc


def wta_product(batch: str="B") -> Construct:
    return Construct([*first_strand(batch),inferred("copied GAAT",revcomp("GAAT")),
                      inferred("copied SMART handle",revcomp(rt.SMART_HANDLE),"tso")],
                     name=f"Drop-seq batch {batch} template-switched cDNA")


def selected_tagmented_strand(batch: str="B",insert_nt: int=34) -> Construct:
    constant="ACGT" if batch=="A" else "AC"
    return Construct([inferred("SMART handle",rt.SMART_HANDLE,"tso"),
        inferred(f"batch {batch} constant",constant),inferred("cell barcode","B"*12,"cbc",placeholder=True),
        inferred("UMI","U"*8,"umi",placeholder=True),inferred("poly(T)","T"*30),
        seg("antisense cDNA","X"*insert_nt,placeholder=True),seg("mosaic end reverse complement",nx.ME_RC,"me"),
        seg("s7 reverse complement",nx.S7_RC,"s7")],name=f"Drop-seq batch {batch} bead-end fragment")


def seq_primers(batch: str="B"):
    return (sp.custom("Read 1",f"Drop-seq batch {batch} custom Read 1",custom_r1(batch),
                      "secondary reconstruction; Drop-seq Table S6 unavailable"),
            sp.NEXTERA["I1"],sp.NEXTERA["R2"])


def final_library(batch: str="B",insert_nt: int=34) -> Construct:
    selected=selected_tagmented_strand(batch,insert_nt)
    lib=Construct([seg("P5",il.P5,"p5"),inferred("spacer",SPACER,"r1"),*list(selected),
                   seg("i7 reverse complement","I"*8,"cbc",placeholder=True),
                   seg("P7 reverse complement",il.P7_RC,"p7")],name=f"Drop-seq batch {batch} library")
    for p in seq_primers(batch):
        if sp.locate(lib,p) is None: raise ValueError(f"{p.name} does not land on {lib.name}")
    return lib


def primer_landings(batch: str="B"):
    lib=final_library(batch); out=[]
    for p in seq_primers(batch):
        h=sp.locate(lib,p)
        if h is None: raise ValueError(f"{p.name} does not land on {lib.name}")
        out.append((p,h))
    return out


def read1_layout(batch: str="B") -> list[tuple[str,str]]:
    return [(x.cycles,x.feature.label) for x in
            sp.feature_spans(final_library(batch),seq_primers(batch),{"Read 1":20})]


def protocol_graph(batch: str="B") -> ProtocolGraph:
    source=("Defining source: https://doi.org/10.1016/j.cell.2015.05.002; "
            "exact oligos reconstructed from the secondary schematic and Seq-Well table")
    graph=ProtocolGraph("Drop-seq",authority=source)
    cells=graph.input("cells","single-cell suspension",MaterialKind.BIOLOGICAL,primary=True)
    beads=graph.input("barcoded_beads","barcoded capture beads",MaterialKind.ADAPTER,
                      state=tuple(bead_primer(batch)),inferred=True,source=source)
    droplets=graph.input("droplet_device","nanoliter droplet generator",MaterialKind.CONDITION)
    captured=graph.product("captured_rna","droplet-confined bead-captured poly(A) RNA",
                           MaterialKind.MOLECULAR,state=capture_scene(batch),inferred=True,
                           source=source)
    empty=graph.product("empty_droplets","unproductive cell/bead loading outcomes",
                        MaterialKind.WASTE)
    graph.react("partition_capture","Co-encapsulate cells and barcoded beads",StepKind.PARTITION,
                ((cells,InputRole.TEMPLATE),(beads,InputRole.STRUCTURAL),
                 (droplets,InputRole.CONDITION)),(captured,empty))
    graph.discard(empty,"droplets without a productive cell-bead pairing do not yield library")
    reverse_transcriptase=graph.input("reverse_transcriptase","reverse transcriptase",
                                     MaterialKind.CATALYST)
    first=graph.product("first_strand",first_strand(batch).name,MaterialKind.MOLECULAR,
                        state=first_strand(batch),inferred=True,source=source)
    graph.react("reverse_transcribe","Synthesize bead-barcoded first-strand cDNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured,InputRole.TEMPLATE),(beads,InputRole.STRUCTURAL),
                 (reverse_transcriptase,InputRole.CATALYST)),(first,))
    tso_input=graph.input("tso","template-switching oligo",MaterialKind.OLIGO,
                          state=tuple(tso()),inferred=True,source=source)
    switched=graph.product("template_switched",wta_product(batch).name,MaterialKind.MOLECULAR,
                           state=wta_product(batch),inferred=True,source=source)
    graph.react("template_switch","Install the distal SMART handle",StepKind.TEMPLATE_SWITCHING,
                ((first,InputRole.TEMPLATE),(tso_input,InputRole.STRUCTURAL),
                 (reverse_transcriptase,InputRole.CATALYST)),(switched,))
    smart=graph.input("smart_primer","SMART amplification primer",MaterialKind.OLIGO,
                      inferred=True,source=source)
    preamp_pol=graph.input("preamp_polymerase","whole-transcriptome PCR polymerase",
                           MaterialKind.CATALYST)
    amplified=graph.product("amplified_cdna","amplified full-length barcoded cDNA",
                            MaterialKind.MOLECULAR,state=wta_product(batch),inferred=True,
                            source=source)
    graph.react("preamp","Amplify template-switched cDNA",StepKind.PCR,
                ((switched,InputRole.TEMPLATE),(smart,InputRole.STRUCTURAL),
                 (preamp_pol,InputRole.CATALYST)),(amplified,),
                geometry=PCRGeometry.symmetric(
                    wta_product(batch),PCRPrimer("SMART amplification primer",rt.SMART_HANDLE),
                    statement="The reconstructed molecule carries the same SMART handle at both amplification boundaries."))
    transposome=graph.input("transposome","Nextera transposome",MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    selected=selected_tagmented_strand(batch)
    tagged=graph.product("selected_tagmented",selected.name,MaterialKind.MOLECULAR,
                         state=selected,inferred=True,source=source)
    waste=graph.product("tagment_waste","tagmentation products without the bead end",
                        MaterialKind.WASTE)
    graph.react("tagment","Tagment cDNA and retain the bead end nearest an s7 insertion",
                StepKind.TAGMENTATION,
                ((amplified,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,waste))
    graph.discard(waste,"fragments lacking the bead-derived handle are not PCR-selected")
    p5=graph.input("p5_primer","P5-SMART selective PCR primer",MaterialKind.OLIGO,
                   state=tuple(p5_hybrid(batch)),inferred=True,source=source)
    p7=graph.input("p7_primer","indexed N7 PCR primer",MaterialKind.OLIGO,
                   inferred=True,source=source)
    polymerase=graph.input("library_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    final=final_library(batch)
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final,
                          inferred=True,source=source)
    p5_sequence=il.P5+SPACER+rt.SMART_HANDLE+("ACGT" if batch=="A" else "AC")
    graph.react("library_pcr","Select the bead-end fragment and add P5/P7",StepKind.PCR,
                ((tagged,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
                 (p7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.exact(
                    selected,PCRPrimer("P5-SMART",p5_sequence,
                                       anneal_nt=len(rt.SMART_HANDLE)+(4 if batch=="A" else 2)),
                    PCRPrimer("indexed N7",il.P7+"N"*8+nx.S7,anneal_nt=len(nx.S7)),
                    statement="The selective PCR geometry is an explicit secondary reconstruction because Drop-seq Table S6 was unavailable."))
    graph.sequence(library,SequencingPlan(
        library=final,primers=seq_primers(batch),acquired_roles=tuple(READ_LENGTHS),
        read_lengths=READ_LENGTHS,
        read_uses={"Read 1":"extract","Index 1 (i7)":"demultiplex","Read 2":"include"}))
    graph.declare_reviewed_scope((
        "cells","barcoded_beads","droplet_device","partition_capture","reverse_transcriptase",
        "reverse_transcribe","tso","template_switch","smart_primer","preamp_polymerase",
        "preamp","transposome","tn5","tagment","p5_primer","p7_primer",
        "library_polymerase","library_pcr","library"),source)
    return graph.require_publishable()
