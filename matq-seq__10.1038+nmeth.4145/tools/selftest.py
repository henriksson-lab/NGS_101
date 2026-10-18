#!/usr/bin/env python3
from pathlib import Path
import sys,re
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import matq_seq as P
from checks import Check,Source,have
from storage import sources_root
c=Check();g=P.protocol_graph()
src=Source(sources_root()/HERE.parent.name/'nmeth.4145_41592_2017_BFnmeth4145_MOESM227_ESM.pdf.txt','published complete MATQ protocol')
if have(c,src,label='published MATQ oligos'):
 text=re.sub(r'\s+','',src.path.read_text())
 for label,seq in (('GAT27dT',P.GAT27_DT),('GAT27-5N3G',P.GAT27_5N3G),('GAT27-5N3T',P.GAT27_5N3T),('GAT21-6N3G',P.GAT21_6N3G),('GAT27',P.GAT27),('3NGAT24',P.THREE_N_GAT24)):
  c(label+' source transcription',seq in text,True)
for step,template in (('wta','duplex'),('conversion','wta_clean'),('library_pcr','ligated')):
 c(step+' uses its actual precursor',g.steps[step].geometry.template is g.materials[template].state,True)
c.report()
