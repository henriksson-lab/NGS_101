#!/usr/bin/env python3
# Base-state transitions, workflow continuity and declared sequencing endpoint are enforced while building.
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from chem import ACE_QC_PRIMERS
from checks import Check
c=Check()
c("publisher Supplementary Protocol 2 QC primer transcription (1502–1505)",
  dict(ACE_QC_PRIMERS), {
    "lambda control forward": "AAAGGCCACTATCAGGCAGCTTTGT",
    "lambda control reverse": "GAAGAAATGTTTCCTGTAATGGAAGATG",
    "T4 control forward": "AGCTTATTTGAAATTAACTGGGTTG",
    "T4 control reverse": "TCTTCTCCAGAATAAATGAATTTACCT",
  })
c.report()
