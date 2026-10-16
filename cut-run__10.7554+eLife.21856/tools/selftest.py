#!/usr/bin/env python3
"""Source readout and optional pellet-fate regressions."""
from pathlib import Path
import re
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from storage import data_root
import cut_run as M
from checks import Check, Source, have
from protocol_graph import Disposition
c=Check()
root=data_root()/"sources"/H.parent.name
source=Source(root/"eLife.21856_v3.xml.txt","tools/get_sources.py --slugs "+H.parent.name)
if have(c,source,label="original eLife v3 methods and appendices"):
    text=source.path.read_text()
    paragraph=text.split("Sequencing libraries were prepared from DNA fragments",1)[1].split("\n",1)[0]
    cycles=int(re.search(r"sequenced for (\d+) cycles in paired-end",paragraph)[1])
    c("both declared acquired reads preserve source paired-cycle count",M.READ_LENGTHS,{"Read 1":cycles,"Read 2":cycles})
    c("source minimum PCR is preserved",M.LIBRARY_MIN_PCR_CYCLES,{"eight":8}[re.search(r"amplifying for (\w+) or more cycles",paragraph)[1]])
    c("pellet QC remains source-positive", "Extract pellet fraction for gel analysis (optional)" in text)
    graph=M.protocol_graph()
    c("optional pellet DNA reaches QC, rather than physical discard",graph.terminals["pellet_qc"].disposition,Disposition.QC)
    c("unassigned source index acquisition is not replaced by a generic eight-cycle run",set(M.READ_LENGTHS),{"Read 1","Read 2"})
c.report()
