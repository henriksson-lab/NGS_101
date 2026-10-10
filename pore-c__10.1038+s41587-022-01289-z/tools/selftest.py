#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
import pore_c as P
from checks import Check
c=Check(); c.section('Pore-C source transcription')
c('NlaIII recognition site',P.NLAIII.site,'CATG')
c('NlaIII produces a 3-prime cohesive end',P.NLAIII.end,'3-prime')
c('example is multiway',sum(s.name.startswith('contact fragment') for s in P.CONCATEMER),3)
c.report()
