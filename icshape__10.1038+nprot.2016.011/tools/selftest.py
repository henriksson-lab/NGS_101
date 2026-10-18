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
source=Source(sources_root()/'icshape__10.1038+nprot.2016.011'/'nprot.2016.011_PMC4896316.xml.txt',"Defining published original; exact manifest identity recorded in source_reviews")
if have(c,source,label="publishedoriginaltranscription"):
    text=source.path.read_text()
    printed=re.search(r'(/5phos/DDDNN[^ ]+)',text).group(1).rstrip(';')
    c('whole interrupted sourceRTprimer',P.RT_ORDER,printed)
    p3=re.search(r'P3-Solexa PCR primer.*?5′-([ACGT ]+)-3′',text).group(1).replace(' ','')
    c('whole sourceP3',P.P3,p3)
    p5=re.search(r'P5-Solexa PCR primer.*?5′-([ACGT ]+)-3′',text,re.S).group(1).replace(' ','')
    c('whole sourceP5',P.P5,p5)
    c('source ligase1 and sequential RT/capture','T4 RNA ligase 1, high concentration' in text and text.index('Resuspend RNA samples from Step 39')<text.index('Add 10 µl of streptavidin bead solution'))
    c('nativeprimer has sourceambiguityclasses and two chemicalspacers',printed.startswith('/5phos/DDDNNAACCNNNN') and printed.count('/iSp18/')==2)
c.report()
