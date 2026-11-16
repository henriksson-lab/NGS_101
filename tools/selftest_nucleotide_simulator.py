#!/usr/bin/env python3
"""Standard-library checks for shared connectivity and explicit fold constraints."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'lib')]
from checks import Check
from lib.nucleotide_simulator import Simulator,Polymerase
from nucleotide_view import NucleotideState

c=Check();s=Simulator();src='synthetic causal control'
ids=s.add('hairpin','ACGTTTTACGT',source=src)
before=s.snapshot();history=len(s.history)
c.raises('mismatched fold rejected',lambda:s.fold('hairpin',0,6,4,source=src))
c('failed fold rolls back every bond and nucleotide',s.snapshot(),before)
c('failed fold leaves no reaction in history',len(s.history),history)
c.raises('zero-length loop rejected',lambda:s.fold('hairpin',0,4,4,source=src))
op=s.fold('hairpin',0,7,4,source=src)
c('fold changes only pairing, never covalent identity',op.before.next,op.after.next)
c('same-strand stem contains four physical base pairs',len(op.after.paired)//2,4)
c.raises('occupied stem cannot be paired twice',lambda:s.fold('hairpin',0,7,4,source=src))
pol=Polymerase('non-displacing DNA polymerase',('DNA',),'DNA')
c('terminal hairpin is not treated as a chemically blocked 3′ end',s.three[ids[-1]],'OH')
bad=NucleotideState(s.snapshot(),(('broken view',(ids[0],ids[2])),))
c.raises('sequence view rejects a fabricated backbone join',bad.rows)
s.unpair(ids,source=src)
c.raises('polymerase cannot extend an unpaired 3′ end',lambda:s.extend('hairpin',pol,source=src))
c.report()
