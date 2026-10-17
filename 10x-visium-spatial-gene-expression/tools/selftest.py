#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import visium
from checks import Check,Source,have,run_common
from storage import data_root
c=Check(); c.section("CG000239 Rev F oligonucleotide appendix")
c("capture poly(dT) has thirty bases",next(s.top for s in visium.FINAL_LIBRARY if s.name=="poly(dT)30"),"T"*30)
c("second-strand and reverse cDNA PCR use the same printed primer",visium.SECOND_STRAND_PRIMER,visium.CDNA_PCR_REVERSE)
c("adapter bottom retains 3-prime T",visium.ADAPTER_BOTTOM[-1],"T")
c("dual 10-cycle sample indices",(len(next(s for s in visium.FINAL_LIBRARY if s.name=="i5")),len(next(s for s in visium.FINAL_LIBRARY if s.name=="i7 reverse complement"))),(10,10))
guide=Source(data_root()/"sources"/H.parent.name/"CG000239_Rev_F.pdf.txt","10x CG000239 Rev F Appendix pp. 69-70")
if have(c,guide,label="CG000239 Rev F slide primer transcription"):
    text=guide.path.read_text(encoding="utf-8")
    c("source prints the modeled slide primer",f"5’-{visium.SLIDE_READ1}-N16-N12-{'T'*30}VN-3’" in text,True)
    c("source prints the modeled cDNA forward primer",f"5’-{visium.CDNA_PCR_FORWARD}-3’" in text,True)
run_common(c); c.report()
