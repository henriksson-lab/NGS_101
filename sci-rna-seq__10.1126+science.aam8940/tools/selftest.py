#!/usr/bin/env python3
"""Check complete sci-RNA-seq templates against the original 2017 preprint."""
from __future__ import annotations
import re
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE)); sys.path.insert(0,str(HERE.parents[1]/"lib"))
from storage import data_root
import sci_rna_seq as S
from checks import Check, Source, have
c=Check(); c.section("Original 2017 preprint oligo templates")
root=data_root()
source=Source(root/"sources"/HERE.parent.name/"104844_v1.full.pdf.txt",
              "Original 2017 Methods; tools/get_sources.py --slugs "+HERE.parent.name)
if have(c,source,label="2017 original preprint"):
    text=source.path.read_text()
    start=text.index("For each well, 1,000-10,000 cells")
    end=text.index("After PCR, samples were pooled",start)
    compact=re.sub(r"\s+","",text[start:end])
    for field in ("[10bpindex]","[i5]","[i7]"):
        compact=compact.replace(field,"N"*10)
    for label,sequence in (("anchored RT primer",S.RT_TEMPLATE),
                           ("complete indexed P5 primer",S.P5_TEMPLATE),
                           ("complete indexed P7 primer",S.P7_TEMPLATE)):
        c(label+" source template",sequence in compact)
c.report()
