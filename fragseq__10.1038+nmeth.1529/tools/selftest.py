#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import chem
from checks import Check
c=Check();g=chem.protocol_graph()
for branch in ('treated','control','pnk'):
    c(branch+' PCR preserves actual first-strand precursor',g.steps[branch+'_pcr'].geometry.template is g.materials[branch+'_cdna'].state,True)
c.report()
