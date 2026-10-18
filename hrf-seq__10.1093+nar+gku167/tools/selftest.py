#!/usr/bin/env python3
from pathlib import Path
import sys,re
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import chem
from checks import Check,Source,have
from storage import sources_root
c=Check();g=chem.protocol_graph()
src=Source(sources_root()/HERE.parent.name/'nar_gku167_supp_supp_gku167_nar-03510-met-f-2013-File009.pdf.txt','published SupplementaryTableS1')
if have(c,src,label='published modified oligos'):
    text=re.sub(r'\s+','',src.path.read_text())
    for name,seq in [('RT_random_primer',chem.HRF_RT_RANDOM),('RT_structure_cassette',chem.HRF_RT_CASSETTE),('LIGATION_ADAPTER_RB','PHO-'+chem.HRF_LIGATION_ADAPTER+'3NHC3'),('PCR_forward',chem.HRF_PCR_FORWARD),*((f'indexed primer {index}',chem._hrf_reverse(index)) for index in chem.HRF_INDEXES)]:
        c(name+' source transcription',seq in text,True)
c('PCR preserves actual ligated single-strand precursor',g.steps['library_pcr'].geometry.template is g.materials['pcr_template'].state,True)
c.report()
