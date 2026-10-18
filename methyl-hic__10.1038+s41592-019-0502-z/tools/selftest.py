#!/usr/bin/env python3
from pathlib import Path
import sys,re
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import methyl_hic as M
from checks import Check,Source,have
from storage import sources_root
c=Check()
s=Source(sources_root()/"methyl-hic__10.1038+s41592-019-0502-z"/"s41592-019-0502-z_41592_2019_502_MOESM1_ESM.pdf.txt","Published MOESM1 bulk sonication/inactivation instructions")
if have(c,s,label="published bulk protocol"):
 text=s.path.read_text()
 enzyme=re.search(r"Inactivate (DpnII) by incubating",text)
 c("published inactivation enzyme",M.DIGEST.enzyme.name,enzyme.group(1))
 cycles=re.search(r"For PE 2x(\d+)",text)
 c("bulk PE read-length reference",M.BULK_READ_NT,int(cycles.group(1)))
c.report()
