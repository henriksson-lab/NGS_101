#!/usr/bin/env python3
"""Source-template transcription checks for sci-RNA-seq."""
from __future__ import annotations
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE)); sys.path.insert(0,str(HERE.parents[1]/"lib"))
import sci_rna_seq as S
from checks import Check
c=Check(); c.section("2017 preprint oligo templates")
c("anchored RT primer template",S.RT_TEMPLATE,"ACGACGCTCTTCCGATCTNNNNNNNNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN")
c("P5 primer template",S.P5_TEMPLATE,"AATGATACGGCGACCACCGAGATCTACACNNNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT")
c("P7 primer template",S.P7_TEMPLATE,"CAAGCAGAAGACGGCATACGAGATNNNNNNNNNNGTCTCGTGGGCTCGG")
c.report()
