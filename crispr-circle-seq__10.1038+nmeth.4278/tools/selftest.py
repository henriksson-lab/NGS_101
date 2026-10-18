#!/usr/bin/env python3
"""Check both modified adapter transcriptions against the published protocol."""
from pathlib import Path
import re
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H), str(H.parents[1]/"lib")]
import crispr_circle as M
from checks import Check, Source, have
from storage import data_root
c=Check()
source=Source(data_root()/"sources"/H.parent.name/"nmeth.4278_41592_2017_BFnmeth4278_MOESM84_ESM.pdf.txt", "published Supplementary Protocol; reacquire via tools/get_sources.py")
if have(c, source, label="published modified adapter sequences"):
    text=source.path.read_text()
    def normalized(s):
        return re.sub(r"\s+", "", s).replace("/ideoxyU/", "U").replace("*", "")
    stem=text.split("CIRCLE-seq Hairpin Adapter",1)[1].split("Annealing Program",1)[0].split("/5Phos/",1)[1]
    library=text.split("* NEBNext Adaptor for Illumina (#E7601A):",1)[1].split("/5Phos/",1)[1].split("-3",1)[0]
    c("oSQT1288 published modified sequence", M.STEM_LOOP, normalized(stem))
    c("E7601A published 64-nt modified reference", M.PUBLISHED_LIBRARY_ADAPTER, normalized(library))
c.report()
