#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import pairedtag as M
from checks import Check,Source,have
from storage import sources_root
c=Check();M.protocol_graph()
source=Source(sources_root()/H.parent.name/"s41592-021-01060-3_supp_NIHMS1660233-supplement-1660233_Supp_Tab1.xlsx.txt","Published PairedTag SupplementaryTable1")
if have(c,source):
    rows={r[0]:r[1] for line in source.path.read_text().split("## sheet 2.")[0].splitlines() if len(r:=line.split("\t"))==2}
    for name,value in (("pMENTs","5Phos/"+M.PMENTS+"ddC"),("Linker-R02",M.LINKER_R02),("Blocker-R02",M.BLOCKER_R02),("Linker-R03",M.LINKER_R03),("Quencher-R03",M.QUENCHER_R03),("Anchor-FokI-GH",M.ANCHOR_WRITTEN[:-1]+"*H"),("P5-FokI",M.P5_FOKI),("P5c-NNDC-FokI","5Phos/"+M.P5C_NNDC),("P5H-FokI",M.P5H_FOKI),("P5Hc-NNDC-FokI","5Phos/"+M.P5HC_NNDC),("PA-F",M.PA_F),("PA-R",M.PA_R),("N5XX",M.N5XX),("P7XX",M.P7XX),("P5 Universal",M.P5_UNIVERSAL[:-1]+"*"+M.P5_UNIVERSAL[-1])):
        c("published "+name,value,rows[name])
    barcode_rows={r[1]:r[2] for line in source.path.read_text().splitlines() if len(r:=line.split("\t"))==3}
    for name,value in (("DNA_#01_RE","/5Phos/"+M.BC1_COMMON+"A"+"CATC"+M.NOTI+M.nx.ME),("RNA_#01_RE","/5Phos/"+M.BC1_COMMON+"T"+"CATC"+M.SBFI+M.RNA_DT),("RNA_#01_NRE","/5Phos/"+M.BC1_COMMON+"T"+"CATC"+M.SBFI+M.RNA_N6)):
        c("published "+name,value,barcode_rows[name])
c.report()
