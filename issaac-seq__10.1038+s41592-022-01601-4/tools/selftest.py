#!/usr/bin/env python3
"""Source-transcription and branch-geometry checks for ISSAAC-seq."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem as P
from checks import Check
import seqprimers

c = Check()
c("ME_S5", P.ME_S5, "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG")
c("ME_S7", P.ME_S7, "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG")
c("ME_bottom bases", P.ME_BOTTOM, "CTGTCTCTTATACACATCT")
c("plate capture handle", P.TRUSEQ_R1_SHORT, "CTACACGACGCTCTTCCGATCT")
c("droplet capture handle", P.TRUSEQ_R2_SHORT, "CAGACGTGTGCTCTTCCGATCT")
c("capture UMI length", P.UMI_LENGTH, 10)
c("capture oligo(dT) length (MOESM4 T30VN)", P.DT_LENGTH, 30)
c("plate run", P.PLATE_READ_LENGTHS,
  {"Read 1": 150, "Index 1 (i7)": 8, "Index 2 (i5)": 8, "Read 2": 150})
c("ATAC primer sites", seqprimers.verify(P.ATAC, P.ATAC_PRIMERS), [])
c("RNA hybrid primer sites", seqprimers.verify(P.RNA, P.RNA_PRIMERS), [])
c("RNA primer families", tuple(p.name.split()[0] for p in P.RNA_PRIMERS),
  ("TruSeq", "Nextera", "TruSeq", "Nextera"))
c.report()
