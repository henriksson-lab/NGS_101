#!/usr/bin/env python3
from pathlib import Path
import csv
import io
import re
import sys
from collections import Counter
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem
from checks import Check,Source,have
from storage import sources_root
c=Check();chem.protocol_graph()
root=sources_root()/H.parent.name
si=Source(root/"pnas.1607794113.sapp.pdf.txt","published SI, Materials and Methods")
if have(c,si,label="published SI modified PCR primers"):
    text="".join(si.path.read_text().split())
    for name,sequence in (("TS-PCR Oligo1",chem.TS_PCR_OLIGO1_WRITTEN),("TS-PCR Oligo2",chem.TS_PCR_OLIGO2_WRITTEN)):
        c(name+" modification transcription",sequence in text,True)
table=Source(root/"pnas.1607794113.sd01.xlsx.txt","published Dataset S1 Table S1")
if have(c,table,label="Dataset S1 library parameters"):
    lines=table.path.read_text().split("## sheet Table S2")[0].splitlines()
    header=next(i for i,line in enumerate(lines) if line.startswith("BotSeqS library ID\t"))
    rows=[r for r in csv.DictReader(io.StringIO("\n".join(lines[header:])),delimiter="\t") if re.fullmatch(r"Bot\d+",r["BotSeqS library ID"])]
    dilutions=Counter(int(r["Final fold-dilution of adapter-ligated molecules"]) for r in rows)
    c("typical dilution from all Table S1 library rows",chem.BOTSEQS_TYPICAL_DILUTION_FOLD,dilutions.most_common(1)[0][0])
    c("source cycle alternatives",chem.PCR_CYCLES,tuple(sorted({int(r["Number of PCR cycles"]) for r in rows})))
    c("source paired reads",{(r["Run type"],r["Number of sequencing cycles"]) for r in rows},{("PE",f"2X{chem.BOTSEQS_READ_BP}")})
c.report()
