#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
from chem import (TRUEPRIME_AMPLIFICATION_C,
                                TRUEPRIME_INACTIVATION_C,
                                TRUEPRIME_PRIMER_LENGTH_NT, workflow)
from checks import Check

c = Check()
# Picher et al., main text/Fig. 3 and Whole genome amplification Methods.
c("principal TthPrimPol primer range", TRUEPRIME_PRIMER_LENGTH_NT, (7, 9))
c("WGA and phi29 inactivation temperatures",
  (TRUEPRIME_AMPLIFICATION_C, TRUEPRIME_INACTIVATION_C), (30, 65))
workflow()
c.report()
