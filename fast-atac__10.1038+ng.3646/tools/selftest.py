#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE),str(HERE.parents[1] / "lib")]
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import chem
from checks import Check,Source,have
from storage import sources_root
import re

c = Check()
c.section("Fast-ATAC source transcription")
c("Corces 2016 Online Methods: 76×8×8×76", chem.READ_LENGTHS,
  {"Read 1": 76, "Index 1 (i7)": 8, "Index 2 (i5)": 8, "Read 2": 76})
c("Both cited dual-index primer families present", tuple(p.role for p in chem.SEQ_PRIMERS),
  ("Read 1", "Index 1 (i7)", "Index 2 (i5)", "Read 2"))
source=Source(sources_root()/"doi_10.1038+nature14590"/
    "nature14590_41586_2015_BFnature14590_MOESM36_ESM.xlsx.txt",
    "Buenrostro2015 Nature14590 published SupplementaryTable1 MOESM36, v2_Ad1/v2_Ad2")
if have(c,source,label="cited dual-index primer families"):
    text=source.path.read_text(encoding="utf-8")
    for family,con in zip(("Ad1","Ad2"),chem.primer_states()):
        members=re.findall(r"^v2_"+family+r"\.\d+_[ACGT]+\t([ACGT]+)",text,re.M)
        pattern="".join(s.top for s in con).replace("N","[ACGT]")
        c("all original v2_"+family+" members match local architecture",
          bool(members) and all(re.fullmatch(pattern,x) is not None for x in members))
chem.protocol_graph()
c.report()
