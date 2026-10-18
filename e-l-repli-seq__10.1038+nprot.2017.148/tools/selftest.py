#!/usr/bin/env python3
from pathlib import Path
import sys,re
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import chem
from checks import Check,Source,have
from storage import sources_root
c=Check();g=chem.protocol_graph()
src=Source(sources_root()/HERE.parent.name/'nprot.2017.148_41596_2018_BFnprot2017148_MOESM34_ESM.pdf.txt','published SupplementaryData1 modified adaptor')
if have(c,src,label='published modified hairpin adaptor'):
    text=re.sub(r'\s+','',src.path.read_text())
    sequence=chem.SUPP_DATA1_ADAPTOR.replace('U','/ideoxyU/')
    c('source adaptor including terminal phosphorothioate','/5Phos/'+sequence[:-1]+'*'+sequence[-1] in text,True)
for prefix in ('early','late'):
    c(prefix+' PCR consumes actual single-stranded IP precursor',g.steps[prefix+'_index_pcr'].geometry.template is g.materials[prefix+'_ip_dna'].state,True)
    c(prefix+' denaturation preserves distinct strand polarities',
      (tuple(g.materials[prefix+'_brdu_strand'].state.strands),tuple(g.materials[prefix+'_unlabelled_strand'].state.strands)),(('nascent',),('partner',)))
    c(prefix+' USER residues enter physical cleanup',prefix+'_user_waste' in {port.material for port in g.steps[prefix+'_ligation_cleanup'].inputs},True)
c.report()
