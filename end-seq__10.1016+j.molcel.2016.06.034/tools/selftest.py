#!/usr/bin/env python3
"""Direct whole modified-hairpin transcription from the defining article."""
from pathlib import Path
import re
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import end_seq as M
from checks import Check, Source, have
from storage import data_root
c=Check()
source=Source(data_root()/"sources"/H.parent.name/"j.molcel.2016.06.034_PMC6299834.xml.txt", "defining article Experimental Procedures END-seq")
if have(c,source,label="published modified END-seq hairpins"):
    text=source.path.read_text()
    for i,adapter in ((1,M.A1),(2,M.A2)):
        match=re.search(r"ENDseq-adaptor-"+str(i)+r", 5′-Phos-(.*?)-3′",text)
        written=re.sub(r"\s+", "", match.group(1)) if match else ""
        tokens=re.findall(r"\[Biotin-dT\]|[ACGTU]",written)
        c("published adaptor"+str(i)+" full bases",adapter.sequence,"".join("T" if x=="[Biotin-dT]" else x for x in tokens))
        c("published adaptor"+str(i)+" biotin-dT positions",adapter.biotin_dt,tuple(j for j,x in enumerate(tokens) if x=="[Biotin-dT]"))
c.report()
