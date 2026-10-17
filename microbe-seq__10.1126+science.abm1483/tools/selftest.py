#!/usr/bin/env python3
from collections import Counter
from pathlib import Path
import re
import sys
H = Path(__file__).resolve().parent; sys.path[:0] = [str(H), str(H.parents[1] / "lib")]
import microbe_seq as M
import storage
from checks import Check, Source, have

c = Check()
c("two-part bead barcode feature", sum(s.feature is not None and s.feature.role == "cell_barcode"
                                       for s in M.FINAL_LIBRARY), 2)
# Source transcription: preprint Table S6 split-and-pool oligo families.
table = Source(storage.sources_root() / "microbe-seq__10.1126+science.abm1483" / "2020.12.14.422699_media-3.xlsx.txt",
               "python3 tools/get_sources.py --slugs microbe-seq__10.1126+science.abm1483 (bioRxiv media-3)")
if have(c, table):
    text = table.path.read_text(encoding="utf-8")
    bc1 = re.findall(rf"Barcode 1 oligo-\d+\t{M.BC1_OLIGO_5P} ([ACGT]+) {M.BC1_OLIGO_3P}", text)
    bc2 = re.findall(rf"Barcode 2 oligo-\d+\t{M.BC2_OLIGO_5P} ([ACGT]+) {M.BC2_OLIGO_3P}", text)
    c("Table S6 barcode-1 oligos per length", dict(sorted(Counter(map(len, bc1)).items())),
      {n: 24 for n in M.BC1_LENGTHS})
    c("Table S6 barcode-2 oligos", (len(bc2), {len(b) for b in bc2}), (384, {M.BC2_LENGTH}))
c.report()
