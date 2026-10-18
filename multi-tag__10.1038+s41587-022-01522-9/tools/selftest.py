#!/usr/bin/env python3
"""MulTI-Tag: build the model and check Supplementary Table 1 transcriptions.

The published MOESM2 workbook is byte-identical (MD5 5cd0f4d7...) to the local
preprint media-2.xlsx, which is the only on-disk copy; tools/source_oligo_check.py
excludes it by filename, so the example rows are checked here against that table.
"""
from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import chem
from checks import Check, Source, have
from storage import sources_root

check = Check()
chem.protocol_graph()

TABLE = Source(sources_root() / HERE.parent.name / "2021.07.08.451691_media-2.xlsx.txt",
               "python3 tools/get_sources.py --slugs " + HERE.parent.name +
               " (published MOESM2 = bioRxiv v2 media-2.xlsx, identical MD5)")
check.section("Supplementary Table 1 example rows")
if have(check, TABLE):
    rows = {}
    for line in TABLE.path.read_text(encoding="utf-8").splitlines():
        cells = line.split("\t")
        if cells and cells[0]:
            rows.setdefault(cells[0], cells)
    check("README prints Tn5MErev as [phos]" + chem.TN5MEREV,
          "Tn5MErev sequence: [phos]" + chem.TN5MEREV in TABLE.path.read_text(encoding="utf-8"))
    for name, seq in (("P5_i5_1", chem.P5_I5_1), ("P5_i5_3",chem.P5_I5_3), ("P5_i5_6",chem.P5_I5_6), ("P7_i7_13", chem.P7_I7_13),
                      ("sci_Ad1.1_TAGATCGC", chem.SCI_AD1_1),
                      ("sci_Ad2.1_TAAGGCGA", chem.SCI_AD2_1)):
        check(f"{name} sequence column", seq in rows.get(name, ()), True)
check.report()
