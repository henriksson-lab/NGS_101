#!/usr/bin/env python3
from pathlib import Path
import sys,re
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import chem
from checks import Check,Source,have
from storage import sources_root
c=Check();g=chem.protocol_graph()
src=Source(sources_root()/HERE.parent.name/'s41592-022-01714-w_PMC9834059.xml.txt','published Nano3P Methods243')
if have(c,src,label='published Nano3P oligos'):
 text=re.sub(r'\s+','',src.path.read_text())
 for label,seq in (('DNA5prime phosphate','/5Phos/'+chem.D_DNA),('RNA3prime spacer',''.join('r'+base for base in chem.R_RNA)+'/3SpC3/'),('CompA',chem.COMPA_DNA)):
  c(label+' source transcription',seq in text,True)
c('RNA templates enter RNase reaction','rt_template_rna' in {port.material for port in g.steps['digest_rna'].inputs},True)
c.report()
