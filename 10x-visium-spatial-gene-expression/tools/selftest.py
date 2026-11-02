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
import visium_rev_a as A
import re
rev_a=Source(sources_root()/H.parent.name/'CG000239_VisiumSpatialGeneExpression_UserGuide_Rev_A.pdf.txt','Exact cited10xCG000239RevA Appendix printedp75')
if have(c,rev_a,label='Rev A predecessor intermediate source transcription'):
    raw=rev_a.path.read_text();text=''.join(raw.split()).replace('’',"'").replace('′',"'")
    patterns=re.findall(r"5'-([ACGT]+)-(N+)-(N+)-(T+)VN-3'",text)
    c('Rev A capture-primer pattern is independently present',bool(patterns),True)
    if patterns:
        handle,cb,umi,tail=patterns[0]
        con=A.full_cdna()
        c('Rev A source handle matches exported construct',con.get('partial Read 1').top,handle)
        c('Rev A source barcode width matches exported feature',len(con.get('spatial barcode')),len(cb))
        c('Rev A source UMI width matches exported feature',len(con.get('UMI')),len(umi))
        c('Rev A source dT width matches exported template',len(con.get('poly(dT)30')),len(tail))
    for key,value in [('forward',A.FORWARD),('reverse',A.REVERSE)]:
        c('Rev A source prints complete cDNA '+key,"5'-"+value+"-3'"in text,True)
    c('Rev A printed amplified TSO-derived end agrees','CCC'+A.revcomp(A.TSO_DNA)in text,True)
    c('Rev A source prints whole modified TSO',"5'-"+A.TSO_DNA+"rGrGrG-3'"in text,True)
    export=A.full_length_export()
    c('Rev A export source stage is before fragmentation',export.material.chemical_state['fragmented'],False)
    c('Rev A export has no newly introduced sample index',export.material.chemical_state['sample_index_present'],False)
c.report()
