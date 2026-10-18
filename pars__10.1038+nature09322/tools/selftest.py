#!/usr/bin/env python3
from pathlib import Path
import sys,re
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import chem
from checks import Check,Source,have
from storage import sources_root
c=Check();g=chem.protocol_graph()
src=Source(sources_root()/HERE.parent.name/'cms_054973.pdf.txt','SOLiD SmallRNAExpressionKit guide4399434 RevC04/2009 AppendixA')
if have(c,src,label='historical kit common primer'):
 text=re.sub(r'\s+','',src.path.read_text())
 c('common P1 source transcription',chem.P1 in text,True)
for prefix in ('v1','s1'):
 c(prefix+' PCR geometry retains actual postRNaseH ss target',g.steps[prefix+'_pcr'].geometry.template is g.materials[prefix+'_pcr_input'].state,True)
c('independent V1/S1 terminal sequencing plans',set(g.sequencing),{'v1_library','s1_library'})
c.report()
