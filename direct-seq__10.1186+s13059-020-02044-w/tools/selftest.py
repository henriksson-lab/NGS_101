#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import direct_seq as P
from checks import Check,Source,have
from storage import sources_root
from chemdraw import Construct,Segment
import nextera
c=Check();P.protocol_graph()
root=sources_root()/H.parent.name
scaffold=Source(root/"s13059-020-02044-w_supp_13059_2020_2044_MOESM3_ESM.xlsx.txt","published Table S2")
if have(c,scaffold,label="published scaffold and tRNA"):
 text=scaffold.path.read_text()
 c("Tail-8A8G capture tract",P.CAPTURE in next(line for line in text.splitlines() if line.startswith("Tail-8A8G\t")),True)
 c("human Gln tRNA",f"Human_tRNA_Gln\t{P.HUMAN_TRNA_GLN}" in text,True)
primers=Source(root/"s13059-020-02044-w_supp_13059_2020_2044_MOESM4_ESM.xlsx.txt","published Table S3d")
if have(c,primers,label="published guide PCR primers"):
 text=primers.path.read_text()
 for name,seq in (("P5_Read1",P.P5_READ1),("tRNA_Read2",P.TRNA_READ2),("P7_Read2_Index1",P.P7_READ2_INDEX1),("P7_Read2_Index2",P.P7_READ2_INDEX2)):
  c("Table S3d "+name,name+"\t"+seq in text,True)
wrong=Construct([Segment(s.name,nextera.READ1_PRIMER,s.tag) if s.name=="TruSeq Read 1 arm" else s for s in P.GUIDE_LIB],"wrong Nextera endpoint")
c.raises("source PCR pair rejects a substituted Nextera Read1 arm",lambda:P.require_guide_primer_sites(wrong))
missing_handle=Construct([s for s in P.FIRST_GUIDE_AMPLICON if s.name!="retained bead Read1 handle"],"incorrectly deleted bead-handle suffix")
c.raises("source primer/bead-handle relation rejects deletion before cell barcode",lambda:P.require_retained_bead_handle(missing_handle))
c.report()
