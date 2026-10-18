#!/usr/bin/env python3
"""Optional checks against explicitlyidentified defining published originals."""
from pathlib import Path
import sys,re
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check,Source,have
from storage import sources_root
c=Check()
source=Source(sources_root()/'paiso-seq__10.1038+s41467-019-13228-9'/'s41467-019-13228-9_supp_41467_2019_13228_MOESM1_ESM.pdf.txt',"Defining published original; exact manifest identity recorded in source_reviews")
if have(c,source,label="publishedoriginaltranscription"):
    text=source.path.read_text()
    table=text.split('Supplementary Table 1.',1)[1].split('Supplementary Table 2.',1)[0]
    compact=re.sub(r'\s+','',table)
    for name,seq in [('TSO',P.TSO_ORDERING),('RTprimer',P.RT_PRIMER),('PCRprimer',P.PCR_PRIMER)]:c(name+' originaltranscription',seq in compact)
    c('whole15modifiedsourceguideoligos',all(P.extension_oligo(i) in compact for i in range(1,16)))
c.report()
