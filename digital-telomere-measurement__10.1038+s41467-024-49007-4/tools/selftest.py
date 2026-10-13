#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import dtm as M
from checks import Check,run_common
c=Check(); c.section("DTM source transcription")
c("seqTether",M.SEQ_TETHER,"AACCTTGGAGATGCACGGAGCAAGCAAT")
c("capture oligo",M.CAPTURE,"TGCTCCGTGCATCTCCAAGGTTCACAAAGACACCGACAACTTTCTTCCTAACC")
c("published capture structure",(len(M.CAPTURE_TETHER_ARM),len(M.CAPTURE_BARCODE),len(M.CAPTURE_TELOMERE_ARM)),(22,25,6))
run_common(c); c.report()
