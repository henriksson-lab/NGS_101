#!/usr/bin/env python3
"""Source ordering and manufacturer-family transcription regressions."""
from pathlib import Path
import os
import re
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as M
from checks import Check, Source, have
c=Check()
root=Path(os.environ.get("CHEM_DATA",H.parents[1]/"_data"))/"sources"/H.parent.name
original=Source(root/"j.cell.2015.03.027_PMC4410947.xml.txt","tools/get_sources.py --slugs "+H.parent.name)
vendor=Source(root/"illumina_TruSeq-SmallRNA.htm.txt","Registered original Illumina TruSeq Small RNA sequence-table topic")
g=M.protocol_graph()
if have(c,original,label="defining 2015 on-bead processing"):
    text=original.path.read_text()
    methods=text.split('mNET-Seq Method and RNA Library Preparation',1)[1].split('Analyses of In',1)[0]
    c("original explicitly identifies phosphatase-minus PNK",'T4 PNK 3′phosphatase minus' in methods)
    c("declared kinase preserves the source enzyme variant",'3-prime phosphatase minus' in g.materials['rna_kinase'].label)
    # Guard the previously reversed physical ordering with an original sequence of
    # clauses, then inspect graph chronology instead of just trusting labels.
    c("source puts bead kinase before extraction",methods.index('T4 PNK')<methods.index('RNA was extracted'))
    order=list(g.steps)
    c("model keeps on-bead kinase before extraction",order.index('phosphorylate')<order.index('extract_rna'))
    c("model phosphorylation substrate is the IP, not extracted RNA",g.steps['phosphorylate'].inputs[0].material,'polii_ip')
    interval=re.search(r'(\d+)–(\d+) nt fragments',methods).groups()
    c("model RNA selection retains source numeric interval",'–'.join(interval)+'-nt' in g.materials['rna_size_rule'].label)
if have(c,vendor,label="whole manufacturer-family adapter sequences"):
    compact=re.sub(r'\s+','',vendor.path.read_text())
    c("whole RNA RA5 manufacturer comparator",M.il.TRUSEQ_SMALL_RNA_RA5.replace('T','U') in compact)
    c("whole DNA RA3 manufacturer comparator",M.il.TRUSEQ_SMALL_RNA_RA3 in compact)
c.report()
