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
s=Source(sources_root()/H.parent.name/'s41592-022-01684-z_41592_2022_1684_MOESM3_ESM.xlsx.txt','Selected published original source; see native migration bounded record')
if have(c,s,label='Original source transcription'):
 compact=re.sub(r'\s+','',s.path.read_text())
 for name,seq in [('sourceprobeAfixedstart',P.PROBE_A_5_SITE+P.PROBE_A_UMI+P.PROBE_A_SPACER),('sourceprobeB3end',P.PROBE_B_SPACER+P.PROBE_B_POLYA+'/3AzideN/'),('sourceU_fwd',P.U_FWD),('sourceP5_family',P.P5_PREFIX+P.P5_INDEX+P.P5_ARM),('sourceP7_N701',P.P7_PREFIX+P.P7_INDEX+P.P7_ARM)]:c(name,seq in compact,True)
c.report()
