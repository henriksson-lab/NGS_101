#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import visium
from checks import Check,Source,have
from storage import sources_root
c=Check()
# Native graph construction enforces continuity and explicit molecular inputs.
visium.protocol_graph()
guide=Source(sources_root()/H.parent.name/"CG000239_Rev_F.pdf.txt","10x CG000239 Rev F Appendix pp. 70–71")
if have(c,guide,label="CG000239 Rev F oligonucleotide transcription"):
    text="".join(guide.path.read_text(encoding="utf-8").split())
    for name,sequence in (("slide primer",visium.SLIDE_READ1+"-N16-N12-"+"T"*30+"VN"),
                          ("cDNA forward",visium.CDNA_PCR_FORWARD),
                          ("cDNA reverse",visium.CDNA_PCR_REVERSE),
                          ("second strand",visium.SECOND_STRAND_PRIMER),
                          ("modified TSO",visium.TSO),
                          ("adapter top",visium.ADAPTER_TOP),
                          ("i5 primer",visium.INDEX_I5_PREFIX+"-N10-"+visium.INDEX_I5_SUFFIX),
                          ("i7 primer",visium.INDEX_I7_PREFIX+"-N10-"+visium.INDEX_I7_SUFFIX)):
        c("source prints "+name,"5’-"+sequence+"-3’" in text,True)
    c("source prints adapter bottom in 3′→5′ orientation", "3’-"+visium.ADAPTER_BOTTOM[::-1]+"-5’" in text,True)
c.report()
