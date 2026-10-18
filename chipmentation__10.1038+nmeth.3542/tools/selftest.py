#!/usr/bin/env python3
from pathlib import Path
import sys,re
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chipmentation as M
from checks import Check,Source,have
from storage import sources_root
c=Check();s=Source(sources_root()/"chipmentation__10.1038+nmeth.3542"/"nmeth.3542_PMC4589892.xml.txt","Published main ChIPmentation and amplification Methods")
if have(c,s):
 text=s.path.read_text()
 match=re.search(r"incubated at (37)°C for (10) minutes",text)
 c("beadbound tagment temperature",M.TAG_TEMP,int(match.group(1)))
 c("beadbound tagment minutes",M.TAG_MINUTES,int(match.group(2)))
 match=re.search(r"(24) cycles of 98°C",text)
 c("analyticalqPCR cyclecount",M.QPCR_CYCLES,int(match.group(1)))
c.report()
