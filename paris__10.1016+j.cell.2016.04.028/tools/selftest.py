#!/usr/bin/env python3
from pathlib import Path
import sys,re
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import chem
from checks import Check,Source,have
from storage import sources_root
c=Check();g=chem.protocol_graph()
src=Source(sources_root()/HERE.parent.name/'PMC4896316_PMC4896316.xml.txt','Flynn2016 original reagent list')
if have(c,src,label='inherited2016 oligos'):
 text=re.sub(r'\s+','',src.path.read_text())
 for name,sequence in [('adapter core',chem.ADAPTER),('P5 Solexa',chem.P5),('P3 Solexa',chem.P3)]:c(name+' source transcription',sequence in text,True)
src=Source(sources_root()/HERE.parent.name/'201804_Protocol_Methods_Mol_Biol.pdf.txt','Lu2018 original comparator oligo table')
if have(c,src,label='2018 comparator P6Tall'):
 from paris_2018_oligos import P6_TALL_V4
 c('P6Tall preserves source leadingCT',P6_TALL_V4 in re.sub(r'\s+','',src.path.read_text()),True)
c('PCR template is the actual purified ss circle',g.steps['library_pcr'].geometry.template is g.materials['clean_circle'].state,True)
c.report()
