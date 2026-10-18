#!/usr/bin/env python3
"""Original published primer-family and deposited-scaffold transcription."""
from pathlib import Path
import sys,re
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/'lib')]
import joung_one_pcr as P
from checks import Check,Source,have
from storage import sources_root,legacy_reference
from plasmid import read_genbank
c=Check();original=Source(sources_root()/H.parent.name/'PMC5526071_PMC5526071.xml.txt','Published manuscript Table3, GeCKO branch')
if have(c,original,label='original Table3'):
 text=original.path.read_text()
 c('all ten forward members',P.FORWARD_PRIMERS,tuple(re.findall(r'^NGS-Lib-Fwd-\d+\t([ACGT]+)\t',text,re.M)))
 c('all eight KO reverse members',P.REVERSE_PRIMERS,tuple(re.findall(r'^NGS-Lib-KO-Rev-\d+\t([ACGT]+)\t',text,re.M)))
map_source=Source(legacy_reference('lenticrispr-v1-screening__10.1126+science.1247005/ref/plasmids/addgene-52961_lentiCRISPRv2.gb'),'Deposited Addgene52961 invariant scaffold component only')
if have(c,map_source,label='invariant deposited scaffold'):
 plasmid=read_genbank(str(map_source.path));f=next(x for x in plasmid.features if x.label=='gRNA scaffold')
 c('reused invariant scaffold matches original feature',P.SCAFFOLD,plasmid.sub(f.start,f.end))
c.report()
