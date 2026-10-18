#!/usr/bin/env python3
"""Construct the protocol-owned declaration and its graph invariants."""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import chem
chem.workflow()

chem.protocol_graph()

from checks import Check, Source, have
from storage import sources_root
import csv
check = Check()
source = Source(sources_root() / HERE.parent.name / "2024.03.17.585235_media-1.xlsx.txt",
                "bioRxiv defining v2 supplement, workbook sheet oligos")
if have(check, source):
    rows = list(csv.reader(source.path.read_text().splitlines(), delimiter="\t"))
    printed = {row[1]: row[2] for row in rows if len(row) > 2 and row[1] in chem.DIRECT_PRIMERS}
    for name, modified in chem.DIRECT_PRIMERS.items():
        check(name + " exact modified sequence in defining workbook", modified, printed.get(name))
check.report()
