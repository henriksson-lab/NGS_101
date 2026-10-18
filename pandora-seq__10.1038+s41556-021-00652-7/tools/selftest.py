#!/usr/bin/env python3
from pathlib import Path
import sys,re
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import chem
from checks import Check,Source,have
from storage import sources_root
c=Check();g=chem.protocol_graph()
c('PCR geometry retains actual first-strand target',g.steps['library_pcr'].geometry.template is g.materials['first_strand'].state,True)
c('RT preserves separate RNA template material',g.materials['rt_rna'].id in g.steps['reverse_transcribe'].outputs,True)
c.report()
