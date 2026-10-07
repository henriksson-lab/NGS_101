"""Construct model for 10x GEM-X Single Cell 3' Gene Expression v4 (CG000731)."""
from __future__ import annotations
import illumina as il
import rt
import seqprimers as sp
from chemdraw import Construct,Scene,Segment

PARTIAL_R1=il.TRUSEQ_READ1[11:]
def seg(name,top,tag=None,**kw): return Segment(name=name,top=top,tag=tag,**kw)
def inferred(name,top,tag=None,**kw): return Segment(name=name,top=top,tag=tag,inferred=True,**kw)
def bead_oligo(): return [seg("Partial Read 1",PARTIAL_R1,"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*12,"umi",placeholder=True),seg("poly(dT)30","T"*30),seg("V","V",placeholder=True),seg("N","N",placeholder=True)]
def inferred_tso(): return [inferred("TSO handle",rt.SMART_HANDLE+"ACAT","tso"),inferred("rGrGrG","GGG","tso")]
def capture_scene():
    mrna=[seg("mRNA","X"*32,placeholder=True),seg("poly(A)","A"*30)]; sc=Scene(); sc.strand("mRNA",mrna,label="mRNA")
    sc.anneal("bead primer",bead_oligo(),to="mRNA",pair=("poly(dT)30","poly(A)"),label="released GEM-X bead primer",unpaired=("V","N")); sc.arrow("bead primer","reverse transcriptase"); return sc
def first_strand(): return Construct([*bead_oligo(),seg("antisense cDNA","X"*34,placeholder=True),inferred("CCC",rt.UNTEMPLATED_TAIL,"tso")],name="v4 first strand")
def switch_scene():
    first=first_strand(); sc=Scene(); sc.strand("first",list(first),label="first strand"); sc.anneal("TSO",inferred_tso(),to="first",pair=("rGrGrG","CCC"),label="TSO",above=True); sc.arrow("first","RT switches template"); return sc
def amplified_cdna(): return Construct([*bead_oligo(),seg("antisense cDNA","X"*34,placeholder=True),inferred("CCC + copied TSO","X"*30,"tso",placeholder=True)],name="v4 amplified cDNA")
def inferred_adapter_top(): return [inferred("Read 2 adapter",il.INDEX1_PRIMER,"r2")]
def inferred_adapter_bottom(): return [inferred("stem complement",il.STEM_COMPLEMENT,"r2"),inferred("3' T","T","r2")]
def adapter_scene():
    sc=Scene(); sc.strand("top",inferred_adapter_top(),label="Read 2 adapter"); sc.anneal("bottom",inferred_adapter_bottom(),to="top",pair=("stem complement","Read 2 adapter"),label="short strand"); sc.mark("bottom","3' T","ligation overhang"); return sc
SEQ_PRIMERS=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["I2"],sp.TRUSEQ["R2"])
def final_library(insert_nt=36):
    lib=Construct([seg("P5",il.P5,"p5"),seg("i5","J"*10,"cbc",placeholder=True),seg("Read 1 arm",il.TRUSEQ_READ1,"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*12,"umi",placeholder=True),seg("poly(dT)30","T"*30),seg("V","V",placeholder=True),seg("N","N",placeholder=True),seg("cDNA insert","X"*insert_nt,placeholder=True),seg("dA junction","A"),seg("Index 1 / Read 2 arm",il.INDEX1_PRIMER,"r2"),seg("i7 reverse complement","I"*10,"cbc",placeholder=True),seg("P7 reverse complement",il.P7_RC,"p7")],name="10x 3' v4 library")
    problems=sp.verify(lib,SEQ_PRIMERS)
    if problems: raise ValueError("invalid v4 library: "+"; ".join(problems))
    return lib
def read1_layout(): return [("1–16","cell barcode"),("17–28","UMI")]
