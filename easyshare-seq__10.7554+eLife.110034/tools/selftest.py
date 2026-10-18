#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from checks import Check,Source,have
from storage import sources_root
c=Check()
import chem as M
M.protocol_graph()
src=Source(sources_root()/H.parent.name/"2024.02.26.581705_media-6.xlsx.txt","Published eLife Supplementary file1, md5-identical cached workbook; MANIFEST.tsv identifies published download")
if have(c,src):
    rows={r[0]:r[1:] for line in src.path.read_text().splitlines() if len(r:=line.split("\t"))>1}
    for name,value in (("Tn5-A",M.TN5_A),("Tn5-reverse","/5Phos/"+M.TN5_REVERSE),("Tn5-B2S","/5Phos/"+M.TN5_B2S),("i7-Tru-Seq-long primer",M.I7_TRUSEQ_LONG),("Nextera N501",M.N501),("Nextera N502",M.N502)):
        c("published workbook "+name,value,rows[name][0])
    for name,value,scaffold in (("Round1_001",M.R1_001,M.R1_001_SCAFFOLD),("Round2_001",M.R2_001,M.R2_001_SCAFFOLD)):
        c("published workbook "+name,(value,scaffold),(rows[name][0].removeprefix("/5Phos/"),rows[name][1]))
    written=rows["Template-Switch oligo"][0]
    c("source TSO RNA residues",M.TSO,written.replace("rU","U").replace("rG","G"))
    c("source RT primer biotin residue and T-tail",M.TN5_B2S+"N"*10+"T"*26+"VN",rows["Reverse-Transcription primer"][0].replace("/5Phos/","").replace("/iBiodT/","T"))
c.report()
