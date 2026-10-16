#!/usr/bin/env python3
"""Dated source architectures, preprint identity and declared kit order."""
from pathlib import Path
import os,re,sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
import chem as M
from checks import Check,Source,have
from protocol_graph import StepKind
c=Check()
b=Path(os.environ.get('CHEM_DATA',H.parents[1]/'_data'))/'sources'/H.parent.name
ion=Source(b/'012UG036V0130_QuantSeq-IonTorrent.pdf.txt','Registered Lexogen June2018 guide')
fwd=Source(b/'015UG009V0260_QuantSeq_Illumina_2021-08-16.pdf.txt','Registered Lexogen August2021 FWD guide')
supp=Source(b/'PMC13309788_supp_gkag652_supplemental_file.pdf.txt','Registered original2026 supplement')
v1=Source(b/'2021.07.14.452379v1.full.pdf.txt','Registered exact bioRxiv version1 original')
if have(c,ion,label='whole original Ion architecture'):
    text=re.sub(r'\s+','',ion.path.read_text())
    for name,seq in (('complete Ion A arm',M.ION_A),('complete Ion P1 arm',M.ION_P1)):
        c(name+' in dated original',seq in text)
if have(c,fwd,label='whole original FWD architecture'):
    text=re.sub(r'\s+','',fwd.path.read_text())
    c('whole FWD source left arm preserved',M.il.P5[:-4]+M.il.TRUSEQ_READ1 in text)
    c('whole FWD opposite adapter preserved',M.revcomp(M.il.TRUSEQ_READ2) in text)
    c('whole source P7-side cap preserved',M.il.P7_RC in text)
    hit=M.sp.locate(M.FWD_LIBRARY,M.FWD_R1,n=100)
    c('variable RNA interval prevents fabricated acquisition into distal tail',hit.reads,'…')
if have(c,supp,label='actual 2026 PCR and readout choices'):
    text=supp.path.read_text()
    methods=text.split('Library preparation was performed by Lexogen GmbH',1)[1].split('Different read depths',1)[0]
    c('PCR count matches original',M.PCR_CYCLES,int(re.search(r'with (\d+) PCR cycles',methods)[1]))
    c('read length matches original',M.READ_CYCLES,int(re.search(r'with (\d+) nt read',methods)[1]))
if have(c,v1,label='actual first preprint identity'):
    text=v1.path.read_text()
    c('version1 is the actual 2021 original, not the 2025 version2','this version posted July 14, 2021' in text[:500])
g=M.protocol_graph();order=list(g.steps)
for prefix in ('ion','fwd'):
    stages=[prefix+'_'+x for x in ('rt','remove_rna','second_strand','purify_ds','pcr','purify_library')]
    positions=[order.index(x) for x in stages]
    c(prefix+' keeps source RNA removal/second strand and both cleanups',positions,sorted(positions))
    # Unknown study revision cannot accidentally become established by a caller
    # omitting an uncertainty argument on a source-reference constructor.
    constructor=M.ion_library if prefix=='ion' else M.fwd_library
    c(prefix+' reference bases are inferred by default',all(s.inferred for s in constructor() if not s.placeholder))
c('platform comparison does not assert a checked physical split',g.steps['platform_variants'].kind,StepKind.MAGIC)
c.report()
