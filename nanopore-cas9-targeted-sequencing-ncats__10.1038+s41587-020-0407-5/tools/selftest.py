#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import ncats as N
from checks import Check
c=Check(); c.section('nCATS source transcription')
c('two 20-nt guide targets are drawn',(len(N.GENOME.get('left Cas9 target')),len(N.GENOME.get('right Cas9 target'))),(20,20))
c('old ends are non-ligatable after CIP',N.SELECTIVE_ENDS.old_ends_phosphorylated,False)
c.report()
