#!/usr/bin/env python3
from pathlib import Path
import sys,re
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import chem
from checks import Check,Source,have
from storage import sources_root
c=Check();g=chem.protocol_graph()
src=Source(sources_root()/HERE.parent.name/'j.cell.2020.08.003_PMC8358727.xml.txt','published2020 KeyResourcesTable')
if have(c,src,label='published SMART oligos'):
 text=re.sub(r'\s+','',src.path.read_text())
 for name,sequence in [('RT primer',chem.RAISIN_RT_PRIMER),('locked TSO',chem.RAISIN_LOCKED_TSO),('PCR primer',chem.RAISIN_PCR_PRIMER)]:c(name+' source transcription',sequence in text,True)
c('Preamp uses actual switched ss target',g.steps['preamp'].geometry.template is g.materials['switched'].state,True)
c('Index PCR uses native tagged substrate',g.steps['library_pcr'].geometry.template is g.materials['tagged'].state,True)
c.report()
