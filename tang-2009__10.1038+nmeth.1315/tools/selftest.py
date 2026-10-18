#!/usr/bin/env python3
"""Optional transcription checks against defining published originals."""
import re
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE));sys.path.insert(0,str(HERE.parents[1]/"lib"))
from checks import Check,Source,have
from storage import sources_root
import tang2009 as Q
check=Check()
source=Source(sources_root()/'tang-2009__10.1038+nmeth.1315'/'nmeth.1315_41592_2009_BFnmeth1315_MOESM306_ESM.xls.txt',"Defining published original; recorded manifest URL and source-review record")
if have(check,source,label="original published oligo transcription"):
    text=source.path.read_text()
    rows={line.split(",",1)[0].strip():line.split(",",1)[1] for line in text.splitlines() if "," in line}
    def bases(name):
        body=rows[name].replace("NH2","")
        sequence="".join(re.findall("[ACGT]+",body))
        return sequence[::-1] if body.lstrip().startswith("3") else sequence
    for label,name,value in (("UP1","Universal Primer 1 (UP1)",Q.UP1),("UP2","Universal Primer 2 (UP2)",Q.UP2),("NH2-UP1","NH2-Universal Primer 1 (NH2-UP1)",Q.AUP1),("NH2-UP2","NH2-Universal Primer 2 (NH2-UP2)",Q.AUP2),("P1upper","P1 Adaptor (Upper strand)",Q.P1_TOP),("P1lower","P1 Adaptor (Lower strand)",Q.P1_BOTTOM),("P2upper","P2 Adaptor (Upper strand)",Q.P2_TOP),("P2lower","P2 Adaptor (Lower strand)",Q.P2_BOTTOM),("PCR1","Library PCR Primer 1",Q.LIBRARY_PCR_1),("PCR2","Library PCR Primer 2",Q.LIBRARY_PCR_2)):
        check(label,value,bases(name))
    check("original NH2 primer/end-bias note retained", "NH2 modification at the 5" in text and "prevents the ligation" in text)
check.report()
