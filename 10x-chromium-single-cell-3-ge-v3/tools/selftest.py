#!/usr/bin/env python3
"""Source-transcription check for the CG000183 gel-bead oligo."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import tenx_v3 as V
from checks import Check
from protocol_graph import SourceDiscrepancy
check=Check()
got="".join(s.top for s in V.bead_oligo()).replace("B","N").replace("U","N")
want="CTACACGACGCTCTTCCGATCT"+"N"*16+"N"*12+"T"*30
check("CG000183 appendix gel-bead oligo transcription", got, want)
check("vendor overview/appendix conflict is attached to the bead primer",
      any(isinstance(note, SourceDiscrepancy) and note.target == "rt_primer"
          for note in V.protocol_graph().annotations), True)
check("CG000183 appendix Capture Seq 1 bead oligo transcription",
      "".join(s.top for s in V.capture_seq1_bead_oligo()).replace("B","N").replace("U","N"),
      "GTCAGATGTGTATAAGAGACAG"+"N"*16+"N"*12+"TTGCTAGGACCGGCCTTAAAGC")
check.report()
