#!/usr/bin/env python3
"""Focused CG000184/CG000197 source-transcription checks."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import tenx_crispr as C
from checks import Check
check=Check(); check.section("10x 3-prime Feature Barcode CRISPR source transcription")
check("Capture Sequence 1",C.CS1,"TTGCTAGGACCGGCCTTAAAGC")
check("Capture Sequence 2",C.CS2,"CCTTAGCCGCTAATAGGTGAGC")
check("Feature cDNA forward",C.FEATURE_CDNA_FORWARD,"GCAGCGTCAGATGTGTATAAGAGACAG")
check("Feature SI reverse",C.FEATURE_SI_REVERSE,"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAAGCAGTGGTATCAACGCAGAG")
check("only variable fields remain placeholders",[s.name for s in C.final_library() if s.placeholder],["cell barcode","UMI","protospacer complement","i7 reverse complement"])
check.report()
