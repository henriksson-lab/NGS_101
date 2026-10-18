#!/usr/bin/env python3
"""Optional checks against explicitlyidentified defining published originals."""
from pathlib import Path
import sys,re
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check,Source,have
from storage import sources_root
c=Check()
source=Source(sources_root()/'cap-seq__10.1371+journal.pgen.1001134'/'PMC2944787_PMC2944787.xml.txt',"Defining published original; exact manifest identity recorded in source_reviews")
if have(c,source,label="publishedoriginaltranscription"):
    text=source.path.read_text()
    washes=re.search(r'washed at (\d+) mM NaCl or (\d+) mM NaCl',text)
    c('source wash concentrations',(P.CAP_HUMAN_WASH_MM,P.CAP_MOUSE_WASH_MM),tuple(map(int,washes.groups())))
    c('source ligation before CAPpurification wording','input DNA was ligated to solexa sequencing adaptors prior to purification' in text)
c.report()
