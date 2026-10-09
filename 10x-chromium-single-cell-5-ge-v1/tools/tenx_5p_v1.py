"""Construct model for 10x Chromium Single Cell 5' Gene Expression v1 (CG000109)."""
from __future__ import annotations
import illumina as il
import rt
import seqprimers as sp
from chemdraw import Construct,Scene,Segment,complement_segments,feature

PARTIAL_R1=il.TRUSEQ_READ1[11:]
SWITCH_SPACER="TTTCTTATAT"
CELL_BARCODE=feature("cell_barcode","cell_barcode","whitelist",whitelist="10x-chromium-5prime-v1")
UMI_FEATURE=feature("umi","umi","random")
SAMPLE_INDEX=feature("sample_index_i7","sample_index","unknown")
FEATURES={"cell barcode":CELL_BARCODE,"UMI":UMI_FEATURE,"i7 reverse complement":SAMPLE_INDEX}
def seg(name,top,tag=None,**kw):
    kw.setdefault("feature",FEATURES.get(name)); return Segment(name=name,top=top,tag=tag,**kw)
def bead_tso(): return [seg("Partial Read 1",PARTIAL_R1,"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*10,"umi",placeholder=True),seg("switch spacer",SWITCH_SPACER,"tso"),seg("rGrGrG","GGG","tso")]
def rt_primer(): return [seg("SMART handle",rt.SMART_HANDLE,"tso"),seg("AC","AC"),seg("poly(dT)30","T"*30),seg("V","V",placeholder=True),seg("N","N",placeholder=True)]
def capture_scene():
    mrna=[seg("mRNA","X"*34,placeholder=True),seg("poly(A)","A"*30)]; sc=Scene(); sc.strand("mRNA",mrna,label="mRNA")
    sc.anneal("RT primer",rt_primer(),to="mRNA",pair=("poly(dT)30","poly(A)"),label="universal RT primer",unpaired=("V","N")); sc.arrow("RT primer","reverse transcriptase"); return sc
def first_strand(): return Construct([*rt_primer(),seg("antisense cDNA","X"*36,placeholder=True),seg("CCC",rt.UNTEMPLATED_TAIL,"tso")],name="v1 first strand")
def switch_scene():
    first=first_strand(); sc=Scene(); sc.strand("first",list(first),label="first strand"); sc.anneal("bead TSO",bead_tso(),to="first",pair=("rGrGrG","CCC"),label="barcoded bead TSO",above=True,mod5="bead"); sc.arrow("first","RT switches template"); return sc
def amplified_cdna(): return Construct([*bead_tso(),seg("5-prime cDNA","X"*40,placeholder=True),*complement_segments(rt_primer())],name="v1 amplified cDNA")
def adapter_top(): return [seg("Read 2 adapter",il.INDEX1_PRIMER,"r2")]
def adapter_bottom(): return [seg("stem complement",il.STEM_COMPLEMENT,"r2"),seg("3' T","T","r2")]
def adapter_scene():
    sc=Scene(); sc.strand("top",adapter_top(),label="Read 2 adapter"); sc.anneal("bottom",adapter_bottom(),to="top",pair=("stem complement","Read 2 adapter"),label="short strand"); sc.mark("bottom","3' T","ligation overhang"); return sc
SEQ_PRIMERS=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["R2"])
def final_library(insert_nt=38):
    lib=Construct([seg("P5",il.P5,"p5"),seg("Read 1 arm",il.TRUSEQ_READ1[4:],"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*10,"umi",placeholder=True),seg("switch spacer",SWITCH_SPACER,"tso"),seg("GGG","GGG","tso"),seg("5-prime cDNA insert","X"*insert_nt,placeholder=True),seg("dA junction","A"),seg("Index 1 / Read 2 arm",il.INDEX1_PRIMER,"r2"),seg("i7 reverse complement","I"*8,"cbc",placeholder=True),seg("P7 reverse complement",il.P7_RC,"p7")],name="10x 5' v1 library")
    problems=sp.verify(lib,SEQ_PRIMERS,required_roles=("Read 1","Index 1 (i7)","Read 2"))
    if problems: raise ValueError("invalid 5' v1 library: "+"; ".join(problems))
    return lib
def read1_layout():
    identifiers=[(x.cycles,x.feature.label) for x in sp.feature_spans(final_library(),SEQ_PRIMERS,{"Read 1":26})]
    return [*identifiers,("27 onward","switch spacer, then transcript 5' end")]
