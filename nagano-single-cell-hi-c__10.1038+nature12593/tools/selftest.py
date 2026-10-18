#!/usr/bin/env python3
from pathlib import Path
import sys,re
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import nagano as M
from checks import Check,Source,have
from storage import sources_root
c=Check()
s=Source(sources_root()/"nagano-single-cell-hi-c__10.1038+nature12593"/"nature12593-s1.pdf.txt","Published original supplement; source ordering strands and terminal modifications")
if have(c,s,label="original supplementary oligos"):
 text=re.sub(r"\s+","",s.path.read_text())
 for name,seq in (("PCR forward",M.PCR_FORWARD),("PCR reverse",M.PCR_REVERSE)):
  c(name+" terminal phosphorothioate transcription",seq[:-1]+"*"+seq[-1] in text,True)
 c("source CAA adapter F",M.il.TRUSEQ_READ1+"CAA*T" in text,True)
 c("source CAA adapter R5P","pTTGAGATCGGAAGAGCGGTTCAGCAGGAATGCCGAG" in text,True)
c.report()
