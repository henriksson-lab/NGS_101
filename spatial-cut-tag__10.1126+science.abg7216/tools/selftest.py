#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import spatial_cut_tag as P
import seqprimers as sp
from checks import Check
c=Check(); c.section("Spatial-CUT&Tag source transcription")
c("representative first barcode oligo", P.se.BARCODE_A_OLIGO_1,
  "/5Phos/AGGCCAGAGCATTCGAACGTGATAGATGTGTATAAGAGACAG")
c("authors' raw Read 1 barcode coordinates",
  [(x.segment, x.cycles) for x in sp.feature_spans(
      P.FINAL_LIBRARY, P.SEQ_PRIMERS,
      {"Read 1": 150, "Index 1 (i7)": 8}) if x.feature.role == "spatial_barcode"],
  [("spatial barcode B", "1–8"), ("spatial barcode A", "39–46")])
c.report()
