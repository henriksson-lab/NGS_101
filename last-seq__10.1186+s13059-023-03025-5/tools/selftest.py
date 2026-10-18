#!/usr/bin/env python3
"""Optional checks against explicitlyidentified defining published originals."""
from pathlib import Path
import sys,re
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import last_seq as P
from checks import Check,Source,have
from storage import sources_root
c=Check()
source=Source(sources_root()/'last-seq__10.1186+s13059-023-03025-5'/'s13059-023-03025-5_supp_13059_2023_3025_MOESM2_ESM.xlsx.txt',"Defining published original; exact manifest identity recorded in source_reviews")
if have(c,source,label="publishedoriginaltranscription"):
    text=source.path.read_text()
    table=text.split('## sheet Table S3',1)[1].split('## sheet Table S4',1)[0]
    rows={line.split('\t')[0]:line.split('\t')[1] for line in table.splitlines() if '\t' in line}
    compact=lambda x:re.sub(r'\s+','',x)
    for name,seq in [('BC2',P.BC2_HAIRPIN),('RandomRT primer6',P.RANDOM_RT),('P5',P.LIB_P5),('P7',P.LIB_P7),('NEBNext Universal',P.NEBNEXT_UNIVERSAL),('NEBNext1',P.NEBNEXT_INDEX_5+'CGTGAT'+P.NEBNEXT_INDEX_3)]:
        c(name,seq,compact(rows[name]).replace('*',''))
    tail=compact(rows['dTrU tail'])
    written='/5Phos/'+P.DTRU_LINKER+P.DTRU_DNA+''.join('r'+x for x in P.DTRU_RNA)+'/3ddC/'
    c('whole physicalDNA/RNA/ddC sourceordering',written,tail)
    c('sourceBC5/7 collision remains explicit',compact(rows['BC5'])==compact(rows['BC7']))
c.report()
