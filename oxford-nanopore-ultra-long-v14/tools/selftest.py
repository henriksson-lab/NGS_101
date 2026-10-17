#!/usr/bin/env python3
from pathlib import Path
import re
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from storage import data_root
import ultra_long as M
from checks import Check, Source, have
c=Check(); c.section("SQK-ULK114 source transcription")
bundle=data_root()/"sources"/H.parent.name
chemistry=Source(bundle/"chemistry-technical-document.html.txt",M.CHEMISTRY_SOURCE)
if have(c,chemistry,label="original adapter table"):
    match=re.search(r"Rapid Adapter \(RA\)\s+Top strand:\s*5[’′']-([ACGT]+)-3",chemistry.path.read_text())
    c("RA top-strand source entry located",match is not None)
    if match:c("published RA strand",M.rapid_adapter_strand().top(),match[1])
protocol=Source(bundle/"ultra-long-dna-sequencing-kit-sqk-ulk114.html.txt",M.SOURCE)
if have(c,protocol,label="original ULK114 protocol"):
    text=protocol.path.read_text()
    dilution=re.search(r"dilute the Fragmentation Mix \(FRA\).*?Fragmentation Mix \(FRA\)\s+(\d+) µl\s+FRA dilution buffer \(FDB\)\s+(\d+) µl",text,re.S|re.I)
    c("FRA dilution source entry located",dilution is not None)
    if dilution:c("FRA:FDB dilution follows source",M.FRA_DILUTION_UL,(int(dilution[1]),int(dilution[2])))
    programme=re.search(r"Room temperature\s+(\d+) minutes\s+(\d+)°C\s+(\d+) minutes\s+On ice\s+Cool on ice for a minimum of (\d+) minutes",text)
    c("tagmentation programme source entry located",programme is not None)
    if programme:
        a,b,d,e=map(int,programme.groups())
        c("tagmentation programme follows source",M.TAGMENTATION_PROGRAM,
          (("room temperature",a),(b,d),("ice",e)))
    attach=re.search(r"Add 5 µl Rapid Adapter \(RA\).*?Incubate the reaction for (\d+) minutes at room temperature",text,re.S)
    c("Rapid Adapter incubation source entry located",attach is not None)
    if attach:c("Rapid Adapter minutes follow source",M.RAPID_ADAPTER_MIN,int(attach[1]))
c.report()
