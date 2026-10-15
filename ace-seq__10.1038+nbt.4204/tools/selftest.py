#!/usr/bin/env python3
# Base-state transitions, workflow continuity and declared sequencing endpoint are enforced while building.
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from chem import ACE_QC_PRIMERS
from checks import Check
c=Check()
c("published ACE-seq control-primer transcription",
  dict(ACE_QC_PRIMERS)["lambda control forward"],
  "AAAGGCCACTATCAGGCAGCTTTGT")
c.report()
