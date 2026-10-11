#!/usr/bin/env python3
from pathlib import Path
import sys
H = Path(__file__).resolve().parent
sys.path[:0] = [str(H), str(H.parents[1] / "lib")]
import chem as P
import seqprimers as sp
from checks import Check

c = Check()
c.section("GeoMx DSP vendor geometry")
c("vendor-published tag regions and lengths",
  [(s.name, s.length_bp) for s in P.dsp.photocleaved_tag()],
  [("proprietary SPR1 PCR handle", 20), ("14-nt UMI", 14),
   ("12-nt RTS target identifier", 12), ("proprietary SPR2 PCR handle", 20)])
c("2x27 identifier cycles",
  [(x.read, x.feature.role, x.cycles) for x in sp.feature_spans(
      P.FINAL_LIBRARY, P.SEQ_PRIMERS, P.READ_LENGTHS)
   if x.read in ("Read 1", "Read 2")],
  [("Read 1", "umi", "1–14"), ("Read 1", "feature_barcode", "15–26"),
   ("Read 2", "feature_barcode", "1–12"), ("Read 2", "umi", "13–26")])
c.report()
