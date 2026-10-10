#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import tell_seq as T
from checks import Check
c=Check(); c.section('TELL-seq source transcription')
c('Index 1 molecular barcode is 18 bases',len(T.FINAL.get('18-nt linked-read barcode')),18)
c('Index 2 sample index is 8 bases',len(T.FINAL.get('8-nt sample index')),8)
c('barcode role is structured',T.FINAL.get('18-nt linked-read barcode').feature.role,'linked_read_barcode')
c.report()
