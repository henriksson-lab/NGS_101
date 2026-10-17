#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import ss_adna as M
from storage import data_root
from checks import Check,Source,have,run_common
c=Check(); run_common(c)
bundle=data_root()/"sources"/"single-stranded-ancient-dna__10.1038+nprot.2013.038"
table=Source(bundle/"nprot.2013.038.pdf.txt","Gansauge & Meyer 2013 corrected article PDF, Table 1")
if have(c,table,label="corrected 2013 Table 1"):
    text=table.path.read_text()
    flat="".join(text.split())
    c("CL9 transcribed",M.CL9 in text)
    c("CL53 with 3' ddC transcribed",M.CL53[:-1]+"-ddC" in text)
    c("CL73 with 5' phosphate and phosphorothioates transcribed",
      "[Phosphate]"+M.CL73[:-4]+"*"+"*".join(M.CL73[-4:]) in text)
    c("corrected CL78 with 5' phosphate transcribed","[Phosphate]"+M.CL78+"[C3Spacer]10[TEG-biotin]" in text)
    c("CL72 transcribed as a standalone row",M.CL72 in text.splitlines())
    c("IS4 transcribed (line-wrapped in Table 1)","IS4P5primerwithoutindex"+M.IS4 in flat)
c.report()
