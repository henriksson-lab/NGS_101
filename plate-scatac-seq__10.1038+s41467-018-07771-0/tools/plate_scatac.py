"""Molecular construct model for Plate_scATAC-seq (Chen et al., 2018)."""
from __future__ import annotations
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct,Scene,Segment,complement_segments,revcomp

def seg(name,top,tag=None,**kw): return Segment(name=name,top=top,tag=tag,**kw)

I7_N701="TCGCCTTA"; I5_S502="CTCTCTAT"
N701=il.P7+I7_N701+nx.S7
S502=il.P5+I5_S502+nx.S5

def n701(): return [seg("P7",il.P7,"p7"),seg("i7",I7_N701,"cbc"),seg("s7",nx.S7,"s7")]
def s502(): return [seg("P5",il.P5,"p5"),seg("i5",I5_S502,"cbc"),seg("s5",nx.S5,"s5")]

def tagmented_scene():
    gap=nx.TAGMENTATION_GAP; core=seg("genomic core","X"*24,placeholder=True)
    top=[seg("s5",nx.S5,"s5"),seg("left ME",nx.ME,"me"),seg("left gap","X"*gap,placeholder=True),core]
    bot=[seg("s7",nx.S7,"s7"),seg("right ME",nx.ME,"me"),seg("right gap","X"*gap,placeholder=True),*complement_segments([core])]
    sc=Scene(); sc.strand("top",top,label="nuclear DNA")
    sc.anneal("bottom",bot,to="top",pair=("genomic core'","genomic core"),label="nuclear DNA")
    sc.anneal("left ME bottom",[seg("left ME'",nx.ME_RC,"me")],to="top",pair=("left ME'","left ME"),label="",mod5="p")
    sc.anneal("right ME bottom",[seg("right ME'",nx.ME_RC,"me")],to="bottom",pair=("right ME'","right ME"),label="",mod5="p",above=True)
    sc.mark("top","left gap","9-nt gap"); sc.mark("bottom","right gap","9-nt gap"); return sc

SEQ_PRIMERS=tuple(sp.NEXTERA[k] for k in ("R1","I1","I2","R2"))
def final_library(insert_nt=34):
    lib=Construct([seg("P5",il.P5,"p5"),seg("i5",I5_S502,"cbc"),seg("s5",nx.S5,"s5"),seg("ME",nx.ME,"me"),
      seg("accessible genomic DNA","X"*insert_nt,placeholder=True),seg("ME reverse complement",nx.ME_RC,"me"),seg("s7 reverse complement",nx.S7_RC,"s7"),
      seg("i7 reverse complement",revcomp(I7_N701),"cbc"),seg("P7 reverse complement",il.P7_RC,"p7")],name="Plate_scATAC-seq library")
    problems=sp.verify(lib,SEQ_PRIMERS)
    if problems: raise ValueError("invalid Plate_scATAC-seq library: "+"; ".join(problems))
    return lib
def primer_landings():
    lib=final_library(); return [(p,sp.locate(lib,p)) for p in SEQ_PRIMERS]
