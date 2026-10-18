#!/usr/bin/env python3
from pathlib import Path
import sys,re
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import chem
from checks import Check,Source,have
from storage import sources_root
c=Check();g=chem.protocol_graph()
c('PCR geometry retains actual purified ss relay aliquot',g.steps['library_pcr'].geometry.template is g.materials['pcr_input'].state,True)
c('ss relay target differs from expected amplified endpoint',g.materials['pcr_input'].state is chem.FINAL_LIBRARY,False)
c.report()
