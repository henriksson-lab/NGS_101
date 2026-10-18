#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from storage import data_root
import zipseq as P
from checks import Check
c=Check();c.section("ZipSeq source transcription");c("zipcode 1",P.ZIP1,"GTAGCAACCACAGATCGCACCCGAGAATTCCATGATGC"+"A"*28);
from checks import Source, have
import re
source=Source(data_root() / 'sources/cite-seq__10.1038+nmeth.4380/nmeth.4380_supp_NIHMS891423-supplement-3.pdf.txt',
              "Original CITE-seq supplemental protocol, 2017-07-03, p.6; fetch via tools/get_sources.py cite-seq")
if have(c,source,label="ZipSeq cited predecessor RPI transcription"):
    source_text=re.sub(r"\s+", "", source.path.read_text())
    c("RPI1 checked against independently extracted original protocol", P.CITE_RPI1 in source_text, True)
c.report()
