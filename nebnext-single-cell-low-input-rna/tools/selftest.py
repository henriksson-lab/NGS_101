#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import low_input_rna as M
from checks import Check,run_common
c=Check(); c.section("E6420 manual oligos"); c("RT primer has the printed anchored 30-T tract",M.RT_PRIMER,"AAGCAGTGGTATCAACGCAGAGTAC"+"T"*30+"V"); c("TSO retains terminal riboguanosines",M.TSO.endswith("rGrGrG"),True); c("the cDNA handles are source-known",all(not s.inferred for s in M.WTA if "handle" in s.name),True); c("all four sequencing reads",len(M.SEQ_PRIMERS),4); run_common(c); c.report()
