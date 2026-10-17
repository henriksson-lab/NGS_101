#!/usr/bin/env python3
"""CoTECH: graph construction plus transcription of the displayed Table S1 oligos."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import chem
import storage
from checks import Check, Source, have

check = Check()
chem.workflow()
chem.protocol_graph()
table = Source(storage.sources_root() / "cotech__10.1038+s41592-021-01129-z"
               / "s41592-021-01129-z_41592_2021_1129_MOESM3_ESM.xlsx.txt",
               "Springer MOESM3 (Supplementary Table 1) of doi:10.1038/s41592-021-01129-z via tools/get_sources.py")
if have(check, table):
    rows = {}
    for line in table.path.read_text(encoding="utf-8").splitlines():
        cells = line.split("\t")
        if len(cells) >= 2:
            rows[cells[0].strip()] = cells[1].strip()
    for name, seq in chem.OLIGOS:
        key = name.split(" (")[0]
        check(f"Table S1 prints {key} as displayed", rows.get(name, rows.get(key)), seq)
check.report()
