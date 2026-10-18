#!/usr/bin/env python3
"""Optional transcription checks against defining published originals."""
import re
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE));sys.path.insert(0,str(HERE.parents[1]/"lib"))
from checks import Check,Source,have
from storage import sources_root
import quartz_seq2 as Q
check=Check()
source=Source(sources_root()/'quartz-seq__10.1186+gb-2013-14-4-r31'/'s13059-018-1407-3_supp_13059_2018_1407_MOESM5_ESM.xlsx.txt',"Defining published original; recorded manifest URL and source-review record")
if have(check,source,label="original published oligo transcription"):
    text=source.path.read_text()
    rows={line.split("\t",1)[0]:line.split("\t")[1] for line in text.splitlines() if "\t" in line}
    for name,value in (("eMDRT0001",Q.RT_PRIMER_1536_FIRST),("Tagging primer",Q.TAGGING_PRIMER),("gM_primer",Q.GM_PRIMER),("rYshapeP5",Q.RYSHAPE_P5),("rYshapeP7LT06",Q.RYSHAPE_P7_LT06),("P5-gMac_hybrid",Q.P5_GMAC),("TPC2",Q.TPC2),("Read1DropQuartz",Q.READ1_DROPQUARTZ)):
        check(name,value,rows[name].replace("*","").upper())
    check("source primer chemical modifications",rows["P5-gMac_hybrid"].count("*")==2 and rows["TPC2"].count("*")==1)
    check("printed v3.2 conflicting headings retained","1536 RT primer with SeqLv15mer" in text and "v3.2 RT primers with 14 mer" in text)
check.report()
