#!/usr/bin/env python3
"""Focused CG000185 source-transcription checks."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import tenx_csp as C
from checks import Check
check=Check(); check.section("CG000185 source transcription")
bead="".join(s.top for s in C.bead_primer()).replace("B","N").replace("U","N")
want="GTCAGATGTGTATAAGAGACAG"+"N"*16+"N"*12+"TTGCTAGGACCGGCCTTAAAGC"
check("gel-bead primer", bead, want)
adt="".join(s.top for s in C.antibody_oligo()).replace("F","N")
want_adt="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"+"N"*10+"N"*15+"N"*9+"GCTTTAAGGCCGGTCCTAGCAA"
check("antibody oligo", adt, want_adt)
check.report()
