#!/usr/bin/env python3
"""Table 1 transcription of placeholder-bearing primers the source-oligo guard cannot check."""
from pathlib import Path
import re
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem
from checks import Check, Source, have
from storage import data_root
c=Check()
text=Source(data_root()/"sources"/H.parent.name/"7651_2021_435.pdf.txt","tools/get_sources.py --slugs "+H.parent.name)
if have(c,text,label="Zheng 2022 Table 1"):
    compact=re.sub(r"\s+","",text.path.read_text())
    c("Table 1 Library PCR Index Primer arms",
      chem.INDEX_5+"[8-bpplatebarcode]"+chem.INDEX_3 in compact)
    c("Table 1 Barcoded RT Primer handle and oligo(dT)VN",
      chem.RT+"[8-bpwellbarcode]"+"N"*10+"T"*25+"VN" in compact)
    c("Table 1 TSO DNA part and riboguanosines",chem.TSO_DNA+"rGrGrG" in compact)
c.report()
