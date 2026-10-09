"""Molecular model for Pi-ATAC-seq; unpublished barcoding-primer regions infer by default."""
from __future__ import annotations
import illumina as il
import nextera as nx
import seqprimers as sp
from chemdraw import Construct,Scene,Segment,complement_segments,feature

def seg(name,top,tag=None,**kw): return Segment(name=name,top=top,tag=tag,**kw)
def inferred(name,top,tag=None,**kw): return Segment(name=name,top=top,tag=tag,inferred=True,**kw)

def tagmented_scene():
 gap=nx.TAGMENTATION_GAP; core=seg("genomic core","X"*24,placeholder=True)
 top=[inferred("s5",nx.S5,"s5"),inferred("left ME",nx.ME,"me"),seg("left gap","X"*gap,placeholder=True),core]
 bot=[inferred("s7",nx.S7,"s7"),inferred("right ME",nx.ME,"me"),seg("right gap","X"*gap,placeholder=True),*complement_segments([core])]
 sc=Scene(); sc.strand("top",top,label="fixed nuclear DNA"); sc.anneal("bottom",bot,to="top",pair=("genomic core'","genomic core"),label="fixed nuclear DNA")
 sc.anneal("left ME bottom",[inferred("left ME'",nx.ME_RC,"me")],to="top",pair=("left ME'","left ME"),label="",mod5="p")
 sc.anneal("right ME bottom",[inferred("right ME'",nx.ME_RC,"me")],to="bottom",pair=("right ME'","right ME"),label="",mod5="p",above=True); return sc

SEQ_PRIMERS=tuple(sp.NEXTERA[k] for k in ("R1","I1","I2","R2"))
def final_library(insert_nt=34):
 lib=Construct([inferred("P5",il.P5,"p5"),inferred("i5","I"*8,"cbc",placeholder=True,feature=feature("sample_i5","sample_index","unknown")),inferred("s5",nx.S5,"s5"),inferred("ME",nx.ME,"me"),
 seg("accessible genomic DNA","X"*insert_nt,placeholder=True),inferred("ME reverse complement",nx.ME_RC,"me"),inferred("s7 reverse complement",nx.S7_RC,"s7"),
 inferred("i7 reverse complement","J"*8,"cbc",placeholder=True,feature=feature("sample_i7","sample_index","unknown")),inferred("P7 reverse complement",il.P7_RC,"p7")],name="inferred Pi-ATAC-seq library")
 problems=sp.verify(lib,SEQ_PRIMERS)
 if problems: raise ValueError("invalid inferred Pi-ATAC-seq library: "+"; ".join(problems))
 return lib
def primer_landings():
 lib=final_library(); return [(p,sp.locate(lib,p)) for p in SEQ_PRIMERS]
