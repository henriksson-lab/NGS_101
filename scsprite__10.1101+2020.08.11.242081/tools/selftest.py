#!/usr/bin/env python3
"""Source checks for reported read minima and the alternate parser prefix."""
from pathlib import Path
import os
import re
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import scsprite as M
from checks import Check,Source,have
c=Check()
bundle=Path(os.environ.get("CHEM_DATA",H.parents[1]/"_data"))/"sources"/"scsprite__10.1101+2020.08.11.242081"
article=Source(bundle/"PMC11588347_PMC11588347.xml.txt",M.SOURCE)
if have(c,article,label="original published read requirements"):
    source=re.search(r"at least (\d+) bp in Read 1.*?and (\d+) bp in Read 2",article.path.read_text())
    c("published read minima located",source is not None)
    if source:c("read minima follow original",tuple(M.READ_MINIMUMS.values()),tuple(map(int,source.groups())))
first=Source(bundle/"author_config_dpm6_y-stag_scSPRITE.txt","author original parser configuration")
second=Source(bundle/"author_config_dpm6_y-stag_scSPRITE2.txt","author alternate parser configuration")
if have(c,first,second,label="original alternate parser prefix"):
    def dpm(p):
        return {fields[1]:fields[2] for line in p.read_text().splitlines() if (fields:=line.split()) and fields[0]=="DPM"}
    a,b=dpm(first.path),dpm(second.path)
    c("all original alternate DPM patterns use the displayed prefix",bool(a) and a.keys()==b.keys() and all(b[k]==M.PARSER_PREFIX+a[k] for k in a))
c.report()
