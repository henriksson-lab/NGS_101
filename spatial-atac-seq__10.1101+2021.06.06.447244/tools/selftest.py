#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem
import seqprimers as sp
from checks import Check
c=Check(); c.section("Spatial-ATAC-seq source transcription (Deng et al. 2022)")
c("Methods Tn5ME-A", chem.TN5_ME_A_WRITTEN, "/5Phos/CATCGGCGTACGACTAGATGTGTATAAGAGACAG")
c("Supplementary Table 3 custom Read 1", chem.CUSTOM_READ_1,
  "GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC")
c("2024 Table 1 Tn5ME-A cannot pair with ligation linker 1",
  chem.FARZAD_ADAPTER_PAIRS_LINKER_1, False)
c("authors' raw Read 1 barcode coordinates",
  [(x.segment, x.cycles) for x in sp.feature_spans(
      chem.FINAL_LIBRARY, chem.SEQ_PRIMERS[:1], {"Read 1": 150})
   if x.feature.role == "spatial_barcode"],
  [("spatial barcode B", "1–8"), ("spatial barcode A", "39–46")])
c.report()
