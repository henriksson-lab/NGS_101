"""Construct model for 10x Chromium Single Cell 3' Gene Expression v3 (CG000183)."""
from __future__ import annotations
import illumina as il
import rt
import seqprimers as sp
from chemdraw import Construct,Scene,Segment,feature

PARTIAL_R1=il.TRUSEQ_READ1[11:]
CELL_BARCODE=feature("cell_barcode","cell_barcode","whitelist",whitelist="10x-chromium-3prime-v3")
UMI_FEATURE=feature("umi","umi","random")
SAMPLE_INDEX=feature("sample_index_i7","sample_index","unknown")
FEATURES={"cell barcode":CELL_BARCODE,"UMI":UMI_FEATURE,"i7 reverse complement":SAMPLE_INDEX}
def seg(name,top,tag=None,**kw):
    kw.setdefault("feature",FEATURES.get(name)); return Segment(name=name,top=top,tag=tag,**kw)
def bead_oligo(): return [seg("Partial Read 1",PARTIAL_R1,"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*12,"umi",placeholder=True),seg("poly(dT)30","T"*30)]
def tso(): return [seg("SMART handle",rt.SMART_HANDLE,"tso"),seg("ACAT","ACAT","tso"),seg("rGrGrG","GGG","tso")]
def capture_scene():
    mrna=[seg("mRNA","X"*32,placeholder=True),seg("poly(A)","A"*30)]; sc=Scene(); sc.strand("mRNA",mrna,label="mRNA")
    sc.anneal("bead primer",bead_oligo(),to="mRNA",pair=("poly(dT)30","poly(A)"),label="released gel-bead primer"); sc.arrow("bead primer","reverse transcriptase"); return sc
def first_strand(): return Construct([*bead_oligo(),seg("antisense cDNA","X"*34,placeholder=True),seg("CCC",rt.UNTEMPLATED_TAIL,"tso")],name="v3 first strand")
def switch_scene():
    first=first_strand(); sc=Scene(); sc.strand("first",list(first),label="first strand"); sc.anneal("TSO",tso(),to="first",pair=("rGrGrG","CCC"),label="TSO",above=True); sc.arrow("first","RT switches template"); return sc
def amplified_cdna(): return Construct([*bead_oligo(),seg("antisense cDNA","X"*34,placeholder=True),seg("CCC",rt.UNTEMPLATED_TAIL,"tso"),seg("copied TSO","X"*27,"tso",placeholder=True)],name="v3 amplified cDNA")
def adapter_top(): return [seg("Partial Read 2",il.INDEX1_PRIMER[:-1],"r2")]
def adapter_bottom(): return [seg("stem complement",il.STEM_COMPLEMENT,"r2"),seg("3' T","T","r2")]
def adapter_scene():
    sc=Scene(); sc.strand("top",adapter_top(),label="Read 2 adapter"); sc.anneal("bottom",adapter_bottom(),to="top",pair=("stem complement","Partial Read 2"),label="short strand"); sc.mark("bottom","3' T","ligation overhang"); return sc
SEQ_PRIMERS=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["R2"])
def final_library(insert_nt=36):
    lib=Construct([seg("P5",il.P5,"p5"),seg("Read 1 arm",il.TRUSEQ_READ1[4:],"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*12,"umi",placeholder=True),seg("poly(dT)30","T"*30),seg("cDNA insert","X"*insert_nt,placeholder=True),seg("dA junction","A"),seg("Index 1 / Read 2 arm",il.INDEX1_PRIMER,"r2"),seg("i7 reverse complement","I"*8,"cbc",placeholder=True),seg("P7 reverse complement",il.P7_RC,"p7")],name="10x 3' v3 library")
    problems=sp.verify(lib,SEQ_PRIMERS,required_roles=("Read 1","Index 1 (i7)","Read 2"))
    if problems: raise ValueError("invalid v3 library: "+"; ".join(problems))
    return lib
def read1_layout(): return [(x.cycles,x.feature.label) for x in sp.feature_spans(final_library(),SEQ_PRIMERS,{"Read 1":28})]
