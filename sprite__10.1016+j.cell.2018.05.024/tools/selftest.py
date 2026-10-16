#!/usr/bin/env python3
"""Transcription guards against the original 2018 article and Table S5."""
from pathlib import Path
import re
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from storage import data_root
import sprite_protocol as M
from checks import Check,Source,have
c=Check()
bundle=data_root()/"sources"/"sprite__10.1016+j.cell.2018.05.024"
paper=Source(bundle/"j.cell.2018.05.024_PMC6548320.xml.txt",M.SOURCE)
if have(c,paper,label="original article PCR primers and read requirements"):
    text=paper.path.read_text()
    paragraph=text.split("Final Library Amplification Primers",1)[1].split("SPRITE Data Processing",1)[0]
    normal=re.sub(r"\s+","",paragraph)
    c("forward PCR primer follows original",M.PCR_FORWARD in normal)
    c("reverse PCR primer retains original index orientation",M.PCR_REVERSE in normal)
    found=re.search(r"at least\s+(\d+)\s*x\s*(\d+)\s*bps",text)
    c("original paired-read minima located",found is not None)
    if found:c("paired-read minima follow original",tuple(M.READ_MINIMUMS.values()),tuple(map(int,found.groups())))
table=Source(bundle/"mmc5.xlsx.txt","author-deposited original 2018 Table S5")
if not table.path.exists():table=Source(bundle/"caltech_mmc5.xlsx.txt",table.origin)
if have(c,table,label="original modified ordering oligos"):
    rows=table.path.read_text().splitlines()
    for label,oligo in M.ORDERING_OLIGOS:
        name=label.split()[-1]
        c(label+" retains original bases and modifications",any(name in r.split("\t") and oligo in r.split("\t") for r in rows))
for filename in ("SPRITE_Protocol_July_2018.pdf.txt", "SPRITE_Protocol.pdf.txt"):
    guide=Source(bundle/filename,"original Guttman lab DNA SPRITE guide (August 9 / April 5, 2018)")
    if have(c,guide,label=filename+" QC primer transcription"):
        text=guide.path.read_text()
        found=re.search(r"DPMQCprimerF\s+DPMQCprimerR\s+5[’′']\s+([ACGT]+)\s+3[’′']\s+5[’′']\s+([ACGT]+)",text)
        c(filename+" named QC entries located",found is not None)
        if found:
            c(filename+" QC primers follow original",(M.DPM_QC_FORWARD,M.DPM_QC_REVERSE),found.groups())
c.report()
