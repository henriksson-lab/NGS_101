#!/usr/bin/env python3
from pathlib import Path
import sys,re
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from checks import Check,Source,have
from storage import sources_root
c=Check()
import chem as M
M.protocol_graph()
source=Source(sources_root()/H.parent.name/"scg_41594_2019_323_MOESM2_ESM.xlsx.txt","Published Table1; Springer original md5-identical perMANIFEST")
if have(c,source):
    rows={r[0]:r[1] for line in source.path.read_text().splitlines() if len(r:=line.split("\t"))==2}
    for name,value in (("pMENTs","5Phos/"+M.MENT_BOTTOM+"ddC"),("Linker-R02",M.ROUND2.linker),("Linker-R03",M.ROUND3.linker),("Linker-R04",M.LINKER_R04),("Blocker-R02",M.ROUND2.blocking_oligo),("Blocker-R03",M.ROUND3.blocking_oligo),("Terminator",M.revcomp(M.LINKER_R04)),("Anchor",M.ANCHOR_WRITTEN),("PA-F",M.PA_F),("PA-R",M.PA_R),("N5XX",M.N5XX),("P7XX",M.P7XX)):
        c("published Table1 "+name,value,rows[name])
barcodes=Source(sources_root()/H.parent.name/"scg_41594_2019_323_MOESM3_ESM.xlsx.txt","Published Table2; Springer original md5-identical perMANIFEST")
if have(c,barcodes):
    rows={r[1]:r[2] for line in barcodes.path.read_text().splitlines() if len(r:=line.split("\t"))==3}
    for name,value in (("DNA_#01_RE","/5Phos/"+M.DNA_BC1_PREFIX+"ATC"+M.NOTI+M.nx.ME),("RNA_#01_RE","/5Phos/"+M.RNA_BC1_PREFIX+"ATC"+M.SBFI+M.RT_DT),("R02_#01","/5Phos/"+M.ROUND2.barcode_5_arm+"AAACCGG"+M.ROUND2.barcode_3_arm),("R03_#01","/5Phos/"+M.ROUND3.barcode_5_arm+"AAACCGG"+M.ROUND3.barcode_3_arm),("R04_#01",M.R4_READ2_ARM+"N"*10+"AAACCGG"+M.R4_LIGATING_ARM)):
        c("published Table2 "+name,value,rows[name])
c.report()
