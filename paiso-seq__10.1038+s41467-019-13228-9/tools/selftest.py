#!/usr/bin/env python3
"""Source-transcription checks for PAIso-seq."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check
c=Check()
c("TSO-derived handle",P.INSERT.segments[0].top,"AAGCAGTGGTATCAACGCAGAGT")
c("the published extension family has fifteen members",len(P.BARCODES),15)
c("bulk BC3 ordering string",P.extension_oligo(3),
  "AAGCAGTGGTATCAACGCAGAGTACTACTAGAGTAGCACTCdUTTTTTTTTdUTTTTTTTT")
c.report()
