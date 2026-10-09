"""Construct model for 10x GEM-X Single Cell 5' Gene Expression v3 (CG000733)."""
from __future__ import annotations
import illumina as il
import rt
import seqprimers as sp
from chemdraw import Construct,Scene,Segment,complement_segments,feature

PARTIAL_R1=il.TRUSEQ_READ1[11:]; SWITCH_SPACER="TTTCTTATAT"
CELL_BARCODE=feature("cell_barcode","cell_barcode","whitelist",whitelist="10x-gem-x-5prime-v3")
UMI_FEATURE=feature("umi","umi","random")
I5_FEATURE=feature("sample_index_i5","sample_index","unknown")
I7_FEATURE=feature("sample_index_i7","sample_index","unknown")
FEATURES={"cell barcode":CELL_BARCODE,"UMI":UMI_FEATURE,"i5":I5_FEATURE,"i7 reverse complement":I7_FEATURE}
def seg(name,top,tag=None,**kw):
    kw.setdefault("feature",FEATURES.get(name)); return Segment(name=name,top=top,tag=tag,**kw)
def inferred(name,top,tag=None,**kw): return Segment(name=name,top=top,tag=tag,inferred=True,**kw)
def bead_tso(): return [seg("Partial Read 1",PARTIAL_R1,"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*12,"umi",placeholder=True),seg("switch spacer",SWITCH_SPACER,"tso"),seg("rGrGrG","GGG","tso")]
def inferred_rt_primer(): return [inferred("RT handle",rt.SMART_HANDLE+"AC","tso"),inferred("poly(dT)30","T"*30),inferred("V","V",placeholder=True),inferred("N","N",placeholder=True)]
def capture_scene():
    mrna=[seg("mRNA","X"*34,placeholder=True),seg("poly(A)","A"*30)]; sc=Scene(); sc.strand("mRNA",mrna,label="mRNA")
    sc.anneal("RT primer",inferred_rt_primer(),to="mRNA",pair=("poly(dT)30","poly(A)"),label="Poly-dT RT Primer B",unpaired=("V","N")); sc.arrow("RT primer","reverse transcriptase"); return sc
def first_strand(): return Construct([*inferred_rt_primer(),seg("antisense cDNA","X"*36,placeholder=True),inferred("CCC",rt.UNTEMPLATED_TAIL,"tso")],name="v3 first strand")
def switch_scene():
    first=first_strand(); sc=Scene(); sc.strand("first",list(first),label="first strand"); sc.anneal("bead TSO",bead_tso(),to="first",pair=("rGrGrG","CCC"),label="barcoded bead TSO",above=True,mod5="bead"); sc.arrow("first","RT switches template"); return sc
def amplified_cdna(): return Construct([*bead_tso(),seg("5-prime cDNA","X"*40,placeholder=True),*complement_segments(inferred_rt_primer())],name="v3 amplified cDNA")
def inferred_adapter_top(): return [inferred("Read 2 adapter",il.INDEX1_PRIMER,"r2")]
def inferred_adapter_bottom(): return [inferred("stem complement",il.STEM_COMPLEMENT,"r2"),inferred("3' T","T","r2")]
def adapter_scene():
    sc=Scene(); sc.strand("top",inferred_adapter_top(),label="Read 2 adapter"); sc.anneal("bottom",inferred_adapter_bottom(),to="top",pair=("stem complement","Read 2 adapter"),label="short strand"); sc.mark("bottom","3' T","ligation overhang"); return sc
SEQ_PRIMERS=(sp.TRUSEQ["R1"],sp.TRUSEQ["I1"],sp.TRUSEQ["I2"],sp.TRUSEQ["R2"])
def final_library(insert_nt=38):
    lib=Construct([seg("P5",il.P5,"p5"),seg("i5","J"*10,"cbc",placeholder=True),seg("Read 1 arm",il.TRUSEQ_READ1,"r1"),seg("cell barcode","B"*16,"cbc",placeholder=True),seg("UMI","U"*12,"umi",placeholder=True),seg("switch spacer",SWITCH_SPACER,"tso"),seg("GGG","GGG","tso"),seg("5-prime cDNA insert","X"*insert_nt,placeholder=True),seg("dA junction","A"),seg("Index 1 / Read 2 arm",il.INDEX1_PRIMER,"r2"),seg("i7 reverse complement","I"*10,"cbc",placeholder=True),seg("P7 reverse complement",il.P7_RC,"p7")],name="10x 5' GEM-X v3 library")
    problems=sp.verify(lib,SEQ_PRIMERS)
    if problems: raise ValueError("invalid 5' v3 library: "+"; ".join(problems))
    return lib
def read1_layout(): return [(x.cycles,x.feature.label) for x in sp.feature_spans(final_library(),SEQ_PRIMERS,{"Read 1":28})]
