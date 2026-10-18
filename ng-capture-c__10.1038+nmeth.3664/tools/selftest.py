#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent
sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from checks import Check
import chem as M
c=Check()
c.section("NG Capture-C published workbook transcription")
c("Hba-1 probe, workbook oligo sheet row94",M.HBA_PROBE,"GATCCATGGTAGCACAGGGCAGCTAAGATGCAAGTCTGAAGGAGGAGTCTGGCGAGCTGCTCCTGCAGTTCCCTGGACCCAGAAGGATGAGCTAGCAGATTCACTTGAGCCAAAGGATTC")
c.report()
