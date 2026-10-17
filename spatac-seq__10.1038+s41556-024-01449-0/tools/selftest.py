#!/usr/bin/env python3
from pathlib import Path
import sys
H = Path(__file__).resolve().parent
sys.path[:0] = [str(H), str(H.parents[1] / "lib")]
import chem as P
import seqprimers as sp
from checks import Check

c = Check()
c.section("SPATAC-seq source transcription")
c("representative modified Tn5 forward oligo", "/5Phos/" + P.spa.TN5_F_D1,
  "/5Phos/AGGCCAGAGCATTCGTATAGCCTAGATGTGTATAAGAGACAG")
c("authors' four cell-barcode read spans",
  [(x.read, x.segment, x.cycles) for x in sp.feature_spans(
      P.FINAL_LIBRARY, P.SEQ_PRIMERS, P.PROCESSING_SPANS)
   if x.feature.role == "cell_barcode"],
  [("Read 1", "cell barcode round 2", "1–8"),
   ("Read 1", "cell barcode round 1 forward", "39–46"),
   ("Read 2", "cell barcode round 3 reverse complement", "1–8"),
   ("Read 2", "cell barcode round 1 reverse complement", "24–31")])
c.report()
