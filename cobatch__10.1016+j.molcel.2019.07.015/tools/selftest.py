#!/usr/bin/env python3
"""Check protocol-owned oligos against the final published Table S1 original."""
from pathlib import Path
import os
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import cobatch as M
from checks import Check, Source, have
c=Check()
bundle=Path(os.environ.get("CHEM_DATA",H.parents[1]/"_data"))/"sources"/"cobatch__10.1016+j.molcel.2019.07.015"
source=Source(bundle/"1-s2.0-S1097276519305453-mmc2.xlsx.txt",M.SOURCE)
if have(c,source,label="original published Table S1 oligos"):
    rows={r.split("\t")[0]:r.split("\t")[1] for r in source.path.read_text().splitlines() if "\t" in r}
    c("complete transferred T5-1 oligo", "".join(s.top for s in M.t5_adapter()),rows["T5-1"])
    c("complete transferred T7-2 oligo", "".join(s.top for s in M.t7_adapter()),rows["T7-2"])
    c("source-modified common ME bottom", "5Phos/"+M.TN5_ME_REV,rows["common annealing primer"])
    c("i5 PCR endpoint arms",M.il.P5+"TAGATCGC"+M.nx.S5,rows["i5-N501"])
    c("i7 PCR endpoint arms",M.il.P7+"TCGCCTTA"+M.nx.S7,rows["i7-N701"])
c.report()
