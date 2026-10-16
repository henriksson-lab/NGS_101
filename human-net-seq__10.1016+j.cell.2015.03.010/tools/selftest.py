#!/usr/bin/env python3
"""Original modified oligos, retained circle junction and bounded read preview."""
from pathlib import Path
import re
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from storage import data_root
import human_net_seq as N
import seqprimers as sp
from checks import Check, Source, have
c=Check()
root=data_root()/"sources"/H.parent.name
original=Source(root/"594009_59c7fec32bcc42cf8c824804675fe319.pdf.txt","Registered author-hosted Mayer2015 Cell PDF; tools/get_sources.py --slugs "+H.parent.name)
predecessor=Source(root/"2012-native-elongating-transcript-sequencing-net-seq.pdf.txt","Registered official Churchman lab 2012 original PDF")
if have(c,original,label="original 2015 modified oligos and readout"):
    text=original.path.read_text()
    methods=text.split("Sequencing Library Preparation for NET-Seq",1)[1].split("Sequencing Library Preparation for Subcellular RNA-Seq",1)[0]
    compact=re.sub(r"\s+","",methods).replace('/iSP18/','/iSp18/')
    c("full preadenylated N6/blocked linker from original",'App'+N.LINKER+'/3ddC' in compact)
    ordered=N.RT_PRIMER_5+'/iSp18/'+N.RT_PRIMER_MID+'/iSp18/'+N.RT_PRIMER_3
    c("entire phosphorylated dual-spacer RT primer from original",'Phos/'+ordered in compact)
    c("full custom read primer from original",N.SEQUENCING in compact)
    c("selected run occurs in original",'NextSeq500(SE75nt)' in compact)
    # The reverse PCR primer binds inside RT5, so circle closure retains its upstream
    # ATC. A generic P7 endpoint previously silently lost that source-known tract.
    rt5=re.search(r'Phos/([ACGT]+)/iSp18/',compact)[1]
    c("final circle-derived tail retains the full original RT5 tract",N.FINAL_LIBRARY.top().endswith(rt5))
if have(c,predecessor,label="inherited original reverse PCR primer"):
    text=predecessor.path.read_text()
    table=text.split('Table 4.14.1 DNA and RNA Oligonucleotide Sequences',1)[1].split('RNA oligonucleotide',1)[0]
    reverse=table.split('oNTI231',1)[1].split('oMHL001',1)[0]
    c("entire reverse primer is in original predecessor table",N.REVERSE_PCR in re.sub(r'\s+','',reverse))
read=sp.locate(N.FINAL_LIBRARY,N.SEQ_PRIMERS[0],n=N.READ_LENGTHS['Read 1']).reads
c("variable biological interval prevents fabricated read-through into RT5",read,'N'*6+'…')
c.report()
