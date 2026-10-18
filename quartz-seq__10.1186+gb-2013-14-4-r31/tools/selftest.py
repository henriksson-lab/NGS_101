#!/usr/bin/env python3
"""Optional transcription checks against defining published originals."""
import re
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE));sys.path.insert(0,str(HERE.parents[1]/"lib"))
from checks import Check,Source,have
from storage import sources_root
import quartz_seq as Q
check=Check()
source=Source(sources_root()/'quartz-seq__10.1186+gb-2013-14-4-r31'/'gb-2013-14-4-r31_PMC4054835.xml.txt',"Defining published original; recorded manifest URL and source-review record")
if have(check,source,label="original published oligo transcription"):
    text=source.path.read_text()
    rows={line.split("\t",1)[0]:line.split("\t",1)[1] for line in text.splitlines() if "\t" in line}
    for name,value in (("RT primer (WTA)",Q.RT_PRIMER),("Tagging primer",Q.TAGGING_PRIMER),("Suppression primer",Q.SUPPRESSION_PRIMER),("TRSU",Q.TRSU),("TRSI-2",Q.TRSI2),("TPC1",Q.TPC1),("TPC2",Q.TPC2)):
        source_bases=rows[name].split(")")[-1].replace("*","").strip().lstrip("-")
        check(name,value,source_bases)
    check("source primer chemical modifications",rows["Suppression primer"].startswith("(NH2)-") and rows["TRSI-2"].startswith("(5'-phosphate)") and "*" in rows["TPC1"] and "*" in rows["TPC2"])
check.report()
