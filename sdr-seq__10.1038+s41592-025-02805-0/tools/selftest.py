#!/usr/bin/env python3
"""Optional transcription against the selected published original."""
from pathlib import Path
import sys,re
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
from checks import Check,Source,have
from storage import sources_root
import chem as P
c=Check()
s=Source(sources_root()/H.parent.name/'s41592-025-02805-0_supp_41592_2025_2805_MOESM3_ESM.xlsx.txt','Selected published original source; see native migration bounded record')
if have(c,s,label='Original source transcription'):
 compact=re.sub(r'\s+','',s.path.read_text())
 for name,seq in [('V2RT_primer1_8',P.RNA_FORWARD+'TCGCCTTA'+'N'*8+'T'*16+'VN'),('Nextera_R1N_S502',P.P5+'CTCTCTAT'+P.R1N_INDEX_ARM),('DNA_N701',P.P7+'TCGCCTTA'+P.R2N_INDEX_ARM),('RNA_N701',P.P7+'TCGCCTTA'+P.R2)]:c(name+' originalV2tabletranscription',seq in compact,True)
c.report()
