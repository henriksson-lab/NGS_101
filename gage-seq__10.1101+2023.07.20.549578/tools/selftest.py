#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from checks import Check,Source,have
from storage import sources_root
c=Check()
import gage_seq as M
M.protocol_graph()
src=Source(sources_root()/H.parent.name/"s41588-024-01745-3_41588_2024_1745_MOESM3_ESM.xlsx.txt","Published2024 Supplementary Table1 MOESM3, NatureGenetics10.1038/s41588-024-01745-3")
if have(c,src):
    rows={r[0].strip():"".join(r[1].split()).upper() for line in src.path.read_text().splitlines() if len(r:=line.split("\t"))>1}
    for name,value in (("RNA-ligation linker-1",M.RNA_LINKER1),("RNA-ligation linker-2",M.RNA_LINKER2_RNA),("Post-TdT-PolyC12-S",M.POST_TDT_POLYC),("Post-TdT-P5(T7)-primer",M.POST_TDT_P5_T7),("Nextera-P5-pre-primer",M.DNA_PRE_F),("Trueseq-P7-pre-P-S",M.DNA_PRE_R),("Nextera-P5-2A",M.DNA_FINAL_F),("Trueseq-8bp-R1",M.DNA_FINAL_R),("Trueseq-F-L-1",M.RNA_FINAL_F),("Nextera-P7-6A",M.RNA_FINAL_R)):
        c("published Table1 "+name,value,rows[name])
    c("source modified dT RT primer",M.RT_DT,rows["BiodT-dt(20)VN RT primer"].replace("/5PHOS/","").replace("/IBIODT/","T"))
    c("source modified N6 RT primer",M.RT_RANDOM,rows["BiodT Radom hexamer RT primer"].replace("/5PHOS/","").replace("/IBIODT/","T"))
c.report()
