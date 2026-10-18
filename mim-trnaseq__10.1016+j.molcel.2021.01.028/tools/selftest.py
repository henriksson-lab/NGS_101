#!/usr/bin/env python3
from pathlib import Path
import sys,re
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import chem
from checks import Check,Source,have
from storage import sources_root
c=Check();g=chem.protocol_graph()
src=Source(sources_root()/HERE.parent.name/'j.molcel.2021.01.028_supp_mmc5.xlsx.txt','published TableS4 modified oligos')
if have(c,src,label='published modified mim oligos'):
 text=re.sub(r'\s+','',src.path.read_text())
 for label,seq in (('I1 with physical ddC','p'+chem.ADAPTER_BASES+'/ddC/'),('RT with RN and nonnucleotide spacer','pRN'+chem.RT_LEFT+'/iSp18/'+chem.RT_BIND),('forward phosphorothioate',chem.MIM_FORWARD[:-1]+'*C'),('reverse phosphorothioate',chem.MIM_REVERSE[:-1]+'*G')):
  c(label+' source transcription',seq in text,True)
c('PCR preserves actual ss-circle aliquot',g.steps['library_pcr'].geometry.template is g.materials['pcr_aliquot'].state,True)
c.report()
