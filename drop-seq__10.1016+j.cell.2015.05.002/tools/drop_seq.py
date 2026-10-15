"""Molecular model for Drop-seq.

Provenance.  The Drop-seq paper (Macosko et al., Cell 2015) prints no oligo: its
Experimental Procedures defer every sequence to Table S6, which sits behind the PMC
download gate.  What the paper does establish 🟢 is the bead itself: beads carrying a
constant PCR handle, a 12-base split-and-pool cell barcode, an 8-base UMI and a T30
capture tail, synthesized by ChemGenes.

The bases below come from the **Seq-Well paper's own Supplementary Table 1**
(Gierahn et al., Nat Methods 14:395, Springer ESM `41592_2017_BFnmeth4179_MOESM237`),
which prints the complete oligo set for the *same* ChemGenes bead, cat. no.
MACOSKO-2011-10, described there as "as described in Macosko et al. Cell 2015":

    1. Barcoded Bead Sequence   5'-Bead-Linker-TTTTTTT + SMART handle + AC + J12 + N8 + T30
    2. Template Switching Oligo SMART handle + GAAT + rGrGrG
    3. SMART PCR Primer         SMART handle
    4. New P5-SMART PCR Hybrid  P5 + GCCTGTCCGCGG + SMART handle + *A*C
    5. Custom Read 1 Primer     GCCTGTCCGCGG + SMART handle + AC

That transfer from the later paper to this one is an **inference**, so every segment it
supplies is created with `inferred=True` and the page draws it dotted.  The `SMART handle`
itself is the canonical `rt.SMART_HANDLE`; P5, P7, s7 and the mosaic end are canonical
`illumina` / `nextera` constants.

The scg_lib_structs drawn page is **not** a source here: it is inspiration only, nothing
is taken from it and agreement with it is not verification.  An earlier revision of this
module drew a second bead batch ("seqA", handle + `ACGT`) and its matching Read 1 primer;
that variant appears in no primary document we hold — only on the drawn page — so it is
no longer modelled.
"""
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
    """Seq-Well Supplementary Table 1 base, carried over to Drop-seq by inference."""
    kw.setdefault("feature", FEATURES.get(name))
    return Segment(name=name, top=top, tag=tag, inferred=True, **kw)


LINKER="T"*7
SPACER="GCCTGTCCGCGG"
BEAD_CONSTANT="AC"      # Seq-Well Supplementary Table 1, row 1


def bead_primer() -> list[Segment]:
    return [inferred("bead linker",LINKER),inferred("SMART handle",rt.SMART_HANDLE,"tso"),
            inferred("bead constant",BEAD_CONSTANT),
            inferred("cell barcode","B"*12,"cbc",placeholder=True),
            inferred("UMI","U"*8,"umi",placeholder=True),inferred("dT30","T"*30)]


def tso() -> list[Segment]:
    return [inferred("SMART handle",rt.SMART_HANDLE,"tso"),inferred("GAAT","GAAT"),
            inferred("rGrGrG","GGG","tso")]


def p5_hybrid() -> list[Segment]:
    return [seg("P5",il.P5,"p5"),inferred("spacer",SPACER,"r1"),
            inferred("SMART handle",rt.SMART_HANDLE,"tso"),
            inferred("bead constant",BEAD_CONSTANT)]


def custom_r1() -> str:
    return SPACER+rt.SMART_HANDLE+BEAD_CONSTANT


def capture_scene() -> Scene:
    mrna=[seg("RNA body","X"*28,placeholder=True),seg("poly(A)","A"*30)]
    sc=Scene(); sc.strand("mRNA",mrna,label="captured mRNA")
    sc.anneal("bead primer",bead_primer(),to="mRNA",pair=("dT30","poly(A)"),
              label="barcoded bead",mod5="bead")
    sc.arrow("bead primer","reverse transcriptase"); return sc


def first_strand() -> Construct:
    return Construct([*bead_primer(),seg("antisense cDNA","X"*28,placeholder=True),
                      inferred("CCC","CCC","tso")],name="Drop-seq first strand")


def template_switch_scene() -> Scene:
    first=first_strand(); sc=Scene(); sc.strand("first strand",list(first),label="first strand",mod5="bead")
    sc.anneal("TSO",tso(),to="first strand",pair=("rGrGrG","CCC"),label="TSO",above=True)
    sc.arrow("first strand","RT switches template"); return sc


def wta_product() -> Construct:
    return Construct([*first_strand(),inferred("copied GAAT",revcomp("GAAT")),
                      inferred("copied SMART handle",revcomp(rt.SMART_HANDLE),"tso")],
                     name="Drop-seq template-switched cDNA")


def selected_tagmented_strand(insert_nt: int=34) -> Construct:
    return Construct([inferred("SMART handle",rt.SMART_HANDLE,"tso"),
        inferred("bead constant",BEAD_CONSTANT),inferred("cell barcode","B"*12,"cbc",placeholder=True),
        inferred("UMI","U"*8,"umi",placeholder=True),inferred("poly(T)","T"*30),
        seg("antisense cDNA","X"*insert_nt,placeholder=True),seg("mosaic end reverse complement",nx.ME_RC,"me"),
        seg("s7 reverse complement",nx.S7_RC,"s7")],name="Drop-seq bead-end fragment")


def seq_primers():
    return (sp.custom("Read 1","Drop-seq custom Read 1",custom_r1(),
                      "published in Seq-Well Supplementary Table 1 for the same ChemGenes "
                      "bead; Drop-seq Table S6 unavailable"),
            sp.NEXTERA["I1"],sp.NEXTERA["R2"])


def final_library(insert_nt: int=34) -> Construct:
    selected=selected_tagmented_strand(insert_nt)
    lib=Construct([seg("P5",il.P5,"p5"),inferred("spacer",SPACER,"r1"),*list(selected),
                   seg("i7 reverse complement","I"*8,"cbc",placeholder=True),
                   seg("P7 reverse complement",il.P7_RC,"p7")],name="Drop-seq library")
    for p in seq_primers():
        if sp.locate(lib,p) is None: raise ValueError(f"{p.name} does not land on {lib.name}")
    return lib


def primer_landings():
    lib=final_library(); out=[]
    for p in seq_primers():
        h=sp.locate(lib,p)
        if h is None: raise ValueError(f"{p.name} does not land on {lib.name}")
        out.append((p,h))
    return out


def read1_layout() -> list[tuple[str,str]]:
    return [(x.cycles,x.feature.label) for x in
            sp.feature_spans(final_library(),seq_primers(),{"Read 1":20})]


def protocol_graph() -> ProtocolGraph:
    source=("Defining source: https://doi.org/10.1016/j.cell.2015.05.002; the exact bases "
            "are those printed in Seq-Well Supplementary Table 1 (doi:10.1038/nmeth.4179) "
            "for the same ChemGenes MACOSKO-2011-10 bead, carried over by inference "
            "because Drop-seq Table S6 is unavailable")
    graph=ProtocolGraph("Drop-seq",authority=source)
    cells=graph.input("cells","single-cell suspension",MaterialKind.BIOLOGICAL,primary=True)
    beads=graph.input("barcoded_beads","barcoded capture beads",MaterialKind.ADAPTER,
                      state=tuple(bead_primer()),inferred=True,source=source)
    droplets=graph.input("droplet_device","nanoliter droplet generator",MaterialKind.CONDITION)
    captured=graph.product("captured_rna","droplet-confined bead-captured poly(A) RNA",
                           MaterialKind.MOLECULAR,state=capture_scene(),inferred=True,
                           source=source)
    empty=graph.product("empty_droplets","unproductive cell/bead loading outcomes",
                        MaterialKind.WASTE)
    graph.react("partition_capture","INFERRED — Co-encapsulate cells and barcoded beads (bead bases from the Seq-Well oligo table)",StepKind.PARTITION,
                ((cells,InputRole.TEMPLATE),(beads,InputRole.STRUCTURAL),
                 (droplets,InputRole.CONDITION)),(captured,empty))
    graph.discard(empty,"droplets without a productive cell-bead pairing do not yield library")
    reverse_transcriptase=graph.input("reverse_transcriptase","reverse transcriptase",
                                     MaterialKind.CATALYST)
    first=graph.product("first_strand",first_strand().name,MaterialKind.MOLECULAR,
                        state=first_strand(),inferred=True,source=source)
    graph.react("reverse_transcribe","INFERRED — Synthesize bead-barcoded first-strand cDNA",
                StepKind.REVERSE_TRANSCRIPTION,
                ((captured,InputRole.TEMPLATE),(beads,InputRole.STRUCTURAL),
                 (reverse_transcriptase,InputRole.CATALYST)),(first,))
    tso_input=graph.input("tso","template-switching oligo",MaterialKind.OLIGO,
                          state=tuple(tso()),inferred=True,source=source)
    switched=graph.product("template_switched",wta_product().name,MaterialKind.MOLECULAR,
                           state=wta_product(),inferred=True,source=source)
    graph.react("template_switch","INFERRED — Install the distal SMART handle",StepKind.TEMPLATE_SWITCHING,
                ((first,InputRole.TEMPLATE),(tso_input,InputRole.STRUCTURAL),
                 (reverse_transcriptase,InputRole.CATALYST)),(switched,))
    smart=graph.input("smart_primer","SMART amplification primer",MaterialKind.OLIGO,
                      inferred=True,source=source)
    preamp_pol=graph.input("preamp_polymerase","whole-transcriptome PCR polymerase",
                           MaterialKind.CATALYST)
    amplified=graph.product("amplified_cdna","amplified full-length barcoded cDNA",
                            MaterialKind.MOLECULAR,state=wta_product(),inferred=True,
                            source=source)
    graph.react("preamp","INFERRED — Amplify template-switched cDNA",StepKind.PCR,
                ((switched,InputRole.TEMPLATE),(smart,InputRole.STRUCTURAL),
                 (preamp_pol,InputRole.CATALYST)),(amplified,),
                geometry=PCRGeometry.symmetric(
                    wta_product(),PCRPrimer("SMART amplification primer",rt.SMART_HANDLE),
                    statement="The molecule carries the same SMART handle at both amplification boundaries, so one primer amplifies it."))
    transposome=graph.input("transposome","Nextera transposome",MaterialKind.ADAPTER)
    tn5=graph.input("tn5","Tn5 transposase",MaterialKind.CATALYST)
    selected=selected_tagmented_strand()
    tagged=graph.product("selected_tagmented",selected.name,MaterialKind.MOLECULAR,
                         state=selected,inferred=True,source=source)
    waste=graph.product("tagment_waste","tagmentation products without the bead end",
                        MaterialKind.WASTE)
    graph.react("tagment","INFERRED — Tagment cDNA and retain the bead end nearest an s7 insertion",
                StepKind.TAGMENTATION,
                ((amplified,InputRole.TEMPLATE),(transposome,InputRole.STRUCTURAL),
                 (tn5,InputRole.CATALYST)),(tagged,waste))
    graph.discard(waste,"fragments lacking the bead-derived handle are not PCR-selected")
    p5=graph.input("p5_primer","P5-SMART selective PCR primer",MaterialKind.OLIGO,
                   state=tuple(p5_hybrid()),inferred=True,source=source)
    p7=graph.input("p7_primer","indexed N7 PCR primer",MaterialKind.OLIGO,
                   inferred=True,source=source)
    polymerase=graph.input("library_polymerase","library PCR polymerase",MaterialKind.CATALYST)
    final=final_library()
    library=graph.product("library",final.name,MaterialKind.LIBRARY,state=final,
                          inferred=True,source=source)
    p5_sequence=il.P5+SPACER+rt.SMART_HANDLE+BEAD_CONSTANT
    graph.react("library_pcr","INFERRED — Select the bead-end fragment and add P5/P7",StepKind.PCR,
                ((tagged,InputRole.TEMPLATE),(p5,InputRole.STRUCTURAL),
                 (p7,InputRole.STRUCTURAL),(polymerase,InputRole.CATALYST)),(library,),
                geometry=PCRGeometry.exact(
                    selected,PCRPrimer("P5-SMART",p5_sequence,
                                       anneal_nt=len(rt.SMART_HANDLE)+len(BEAD_CONSTANT)),
                    PCRPrimer("indexed N7",il.P7+"N"*8+nx.S7,anneal_nt=len(nx.S7)),
                    statement="The two PCR primers are the ones Seq-Well Supplementary Table 1 prints for this bead; their transfer to Drop-seq is inferred because Drop-seq Table S6 is unavailable."))
    graph.sequence(library,SequencingPlan(
        library=final,primers=seq_primers(),acquired_roles=tuple(READ_LENGTHS),
        read_lengths=READ_LENGTHS,
        read_uses={"Read 1":"extract","Index 1 (i7)":"demultiplex","Read 2":"include"}))
    graph.declare_reviewed_scope((
        "cells","barcoded_beads","droplet_device","partition_capture","reverse_transcriptase",
        "reverse_transcribe","tso","template_switch","smart_primer","preamp_polymerase",
        "preamp","transposome","tn5","tagment","p5_primer","p7_primer",
        "library_polymerase","library_pcr","library"),source)
    return graph.require_publishable()
