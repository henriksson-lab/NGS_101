#!/usr/bin/env python3
from pathlib import Path
import re
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from storage import data_root
import rapid_v14 as M
from checks import Check, Source, have
c=Check(); c.section("SQK-RAD114 source transcription")
bundle=data_root()/"sources"/H.parent.name
chemistry=Source(bundle/"chemistry-technical-document.html.txt",M.CHEMISTRY_SOURCE)
if have(c,chemistry,label="original adapter table"):
    match=re.search(r"Rapid Adapter \(RA\)\s+Top strand:\s*5[’′']-([ACGT]+)-3",chemistry.path.read_text())
    c("RA top-strand source entry located",match is not None)
    if match:c("published RA strand",M.rapid_adapter_strand().top(),match[1])
protocol=Source(bundle/"rapid-sequencing-sqk-rad114.html.txt",M.SOURCE)
if have(c,protocol,label="original Rapid V14 protocol"):
    text=protocol.path.read_text()
    programme=re.search(r"Incubate the tube at (\d+)°C for (\d+) minutes and then at (\d+)°C for (\d+) minutes",text)
    c("fragmentation programme source entry located",programme is not None)
    if programme:
        a,b,d,e=map(int,programme.groups())
        c("fragmentation programme follows source",M.FRAGMENTATION_PROGRAM,((a,b),(d,e)))
    attachment=re.search(r"Incubate the reaction for (\d+) minutes at room temperature",text)
    c("attachment source entry located",attachment is not None)
    if attachment:c("Rapid Adapter attachment minutes follow source",M.RAPID_ATTACHMENT_MIN,int(attachment[1]))
    dilution=re.search(r"dilute the Rapid Adapter \(RA\).*?Rapid Adapter \(RA\)\s+([\d.]+) [µμ]l\s+Adapter Buffer \(ADB\)\s+([\d.]+) [µμ]l",text,re.S)
    c("RA dilution source entry located",dilution is not None)
    if dilution:c("RA:ADB dilution follows source",M.RA_DILUTION_UL,(float(dilution[1]),float(dilution[2])))
c.report()
