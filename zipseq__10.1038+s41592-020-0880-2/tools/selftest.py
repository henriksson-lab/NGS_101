#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from storage import data_root
import zipseq as P
from checks import Check
c=Check();c.section("ZipSeq source transcription");c("zipcode 1",P.ZIP1,"GTAGCAACCACAGATCGCACCCGAGAATTCCATGATGC"+"A"*28);
from chemdraw import Construct, Segment, revcomp
from checks import Source, have
import illumina, nextera, re
wrong_read1 = Construct([
    Segment(s.name, nextera.READ1_PRIMER, s.tag) if s.name == "TruSeq Read 1 arm" else s
    for s in P.ZIP_LIB], "wrong Nextera zipcode endpoint")
c.raises("cited P5 primer rejects a Nextera barcode-end substitution",
         lambda: P.require_zipcode_primer_sites(wrong_read1))
wrong_read2 = Construct([
    Segment(s.name, revcomp(illumina.TRUSEQ_READ2), s.tag)
    if s.name == "RPI-derived extension of Small RNA arm" else s
    for s in P.ZIP_LIB], "wrong ordinary TruSeq zipcode endpoint")
c.raises("cited RPI architecture rejects ordinary TruSeq Read2 substitution",
         lambda: P.require_zipcode_primer_sites(wrong_read2))
source=Source(data_root() / 'sources/cite-seq__10.1038+nmeth.4380/nmeth.4380_supp_NIHMS891423-supplement-3.pdf.txt',
              "Original CITE-seq supplemental protocol, 2017-07-03, p.6; fetch via tools/get_sources.py cite-seq")
if have(c,source,label="ZipSeq cited predecessor RPI transcription"):
    source_text=re.sub(r"\s+", "", source.path.read_text())
    c("RPI1 checked against independently extracted original protocol", P.CITE_RPI1 in source_text, True)
c.report()
