#!/usr/bin/env python3
"""Focused CG000185 source-transcription checks."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import tenx_csp as C
bead="".join(s.top for s in C.bead_primer()).replace("B","N").replace("U","N")
want="GTCAGATGTGTATAAGAGACAG"+"N"*16+"N"*12+"TTGCTAGGACCGGCCTTAAAGC"
if bead != want: raise SystemExit("FAIL: CG000185 gel-bead primer transcription")
adt="".join(s.top for s in C.antibody_oligo()).replace("F","N")
want_adt="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"+"N"*10+"N"*15+"N"*9+"GCTTTAAGGCCGGTCCTAGCAA"
if adt != want_adt: raise SystemExit("FAIL: CG000185 antibody oligo transcription")
raise SystemExit(0)
