#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
from chem import (BOTSEQS_INSERT_BP, BOTSEQS_READ_BP,
                                BOTSEQS_TYPICAL_DILUTION_FOLD, workflow)
from checks import Check

c = Check()
# Hoang et al., BotSeqS library preparation and Results.
c("typical dilution after adapter ligation", BOTSEQS_TYPICAL_DILUTION_FOLD, 100_000)
c("selected insert and paired-read lengths", (BOTSEQS_INSERT_BP, BOTSEQS_READ_BP),
  (350, 100))
workflow()
c.report()
