#!/usr/bin/env python3
"""Source-transcription check for the representative classic BLISS adapter A1.

Construction is owned locally; exact finishing geometry remains unresolved.
"""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from storage import data_root
import bliss_chem as B
from checks import Check, Source, have

check = Check()
data1 = Source(data_root() / "sources" / HERE.parent.name / "ncomms15058_supp_ncomms15058-s2.xlsx.txt",
               "published Supplementary Data 1; reacquire via tools/get_sources.py")
if have(check, data1, label="Supplementary Data 1 A1 transcription"):
    rows = [line.split("\t") for line in data1.path.read_text().splitlines()]
    a1 = next(i for i, r in enumerate(rows) if r[0] == "A1" and r[1] == "Bottom")
    check("Supplementary Data 1 classic A1 bottom oligo", B.A1_BOTTOM, rows[a1][2])
    check("Supplementary Data 1 classic A1 top oligo", B.A1_TOP, rows[a1 + 1][2])
check.report()
