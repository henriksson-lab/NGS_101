#!/usr/bin/env python3
"""Cell 2015 Table S6 transcription for oligos the graph does not carry as exact states.

Pairing and primer-site guarantees live in the constructors used by build_page.py; the
graph's exact oligo states are checked by tools/selftest_source_oligos.py.
"""
from pathlib import Path
import re
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
from storage import data_root
from checks import Check,Source,have
import drop_seq as D

check=Check()
s=Source(data_root()/"sources"/HERE.parent.name/"j.cell.2015.05.002_1-s2.0-S0092867415005498-mmc1.pdf.txt",
         "Macosko et al. 2015 Cell Supplemental Information, "
         "https://ars.els-cdn.com/content/image/1-s2.0-S0092867415005498-mmc1.pdf (Table S6)")
if have(check,s,label="Cell 2015 Table S6"):
    text=re.sub(r"\s+","",s.path.read_text(encoding="utf-8"))
    for v in ("SeqA","SeqB"):
        check(f"Read1CustomSeq{v[-1]} printed",D.custom_r1(v) in text,True)
        bead="".join(x.top if not x.placeholder else {"cell barcode":"J","UMI":"N"}[x.name]*len(x.top)
                     for x in D.bead_primer(v))
        check(f"Barcoded Bead {v} printed","Linker-"+bead+"-3" in text,True)
check.report()
