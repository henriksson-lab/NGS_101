#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import native_barcoding as M
from checks import Check
c=Check(); c.section("SQK-NBD114.24 source transcription")
c("NB01", M.NB01, "CACAAAGACACCGACAACTTTCTT")
c("CHTD barcode flanks, forward", M.NB_FLANK_FORWARD, ("AAGGTTAA", "CAGCACCT"))
c("CHTD barcode flanks, reverse", M.NB_FLANK_REVERSE, ("GGTGCTG", "TTAACCTTAGCAAT"))
c("published barcode count", len(M.NATIVE_BARCODES), 24)
c("end-prep programme", M.END_PREP_PROGRAM, ((20, 5), (65, 5)))
c("barcode ligation minutes", M.BARCODE_LIGATION_MIN, 20)
c("Native Adapter ligation minutes", M.NATIVE_ADAPTER_LIGATION_MIN, 20)
c.report()
