#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import zipseq as P
from checks import Check
c=Check();c.section("ZipSeq source transcription");c("zipcode 1",P.ZIP1,"GTAGCAACCACAGATCGCACCCGAGAATTCCATGATGCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");c.report()
