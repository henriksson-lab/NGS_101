#!/usr/bin/env python3
"""Source-transcription checks for scRICA-seq."""
import re
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from storage import data_root
import chem as P
from checks import Check, Source, have
c=Check()
c("source-distinct TSS/TES size fractions and accessibility output",
  tuple(row[0] for row in P.FINAL_LIBRARIES),
  ("TSS long-fragment library","TSS short-fragment library",
   "TES long-fragment library","TES short-fragment library","Accessibility library"))
c("linker is the two exact end homologies", P.LINKER,
  "ACTCTGCGTTGATACCACTGCTTGCAATGAAGTCGCAGGGTTG")
c("source-printed TSS enrichment primer", P.I7_ANTI_TAG1,
  "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGCGCGTGCCGTATAGCGTTAACT")
c("source-printed terminal primer pools", (P.INDEX_N7, P.INDEX_N5),
  ("CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG",
   "AATGATACGGCGACCACCGAGATCTACACNNNNNNNNTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"))
source_dir=data_root()/"sources"/H.parent.name
oligo_table=Source(source_dir/"s41467-024-52335-0_supp_41467_2024_52335_MOESM5_ESM.xlsx.txt",
                   "Jiang et al. 2024 Supplementary Data 2; fetch with tools/get_sources.py")
if have(c, oligo_table, label="scRICA-seq low-throughput oligo table"):
    rows=oligo_table.path.read_text(encoding="utf-8").splitlines()
    oligo_row=next((line for line in rows if line.startswith("Oligo(dT) primer\t")),"")
    match=re.search(r"\(N\)(\d+)\(N\)(\d+)",oligo_row)
    c("printed low-throughput barcode and UMI lengths",
      tuple(map(int,match.groups())) if match else None,
      (len(P.CIRCLE_TEMPLATE.get("cell barcode")),len(P.CIRCLE_TEMPLATE.get("UMI"))))
c("low-throughput ATAC branch has no molecular cell barcode",
  any(seg.feature and seg.feature.role=="cell_barcode" for seg in P.ATAC.segments),False)
c.report()
