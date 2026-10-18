#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
from checks import Check,Source,have
from storage import sources_root
c=Check()
base=sources_root()/'sureplex'
s=Source(base/'illumina_sureplex-protocol-reference-guide-15053626-02.pdf.txt','SurePlex15053626v02 original')
if have(c,s):
 t=''.join(s.path.read_text().split());c('source SurePlex identity','15053626v02' in t and 'September2020' in t,True)
s=Source(base/'illumina_veriseq-pgs-library-prep-reference-guide-15052877-05.pdf.txt','VeriSeq15052877v05 original')
if have(c,s):
 t=''.join(s.path.read_text().split());c('source VeriSeq identity','15052877v05' in t and 'August2021' in t,True)
 c('native normalized library is source ssDNA','normalizedlibrariesconsistofsingle-strandedDNA(ssDNA)' in t,True)
s=Source(base/'illumina_veriseq-pgs-guidance-note-1000000138837-00.pdf.txt','GuidanceNote1000000138837v00 original')
if have(c,s):
 t=''.join(s.path.read_text().split());c('control unmapped reads are computationally filtered','donotmapagainstthehumangenomeandarefilteredoutbythesystem' in t,True)
 c('source36-cycle indexed short-read recipe','36cyclesingle-anddual-indexreadsequencingrunonly' in t,True)
c.report()
