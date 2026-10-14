#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import telomere_profiling as M
from checks import Check,run_common
c=Check(); c.section("Telomere Profiling source transcription")
c("NB50 complete TeloTag",M.TELOTAG_TOPS["NB50"],"CCCTCCGATAATGGACTTTGGTAACTTCCTGCGTTGCGTACAGCAATCAGGCGCGCCGAATTCCGATCGTTAATTAATGGTGGACTCTATGACCGTTCA")
c("NB50.C1 complete splint",M.splint_family("NB50")[0],"GAGTCCACCATTAATTAACGATCGGAATTCGGCGCGCCTGATTGCTGTACGCAACGCAGGAAGTTACCAAAGTCCATTATCGGAGGGACCCTAACCCTAACCCTA")
c("capture has release step",len(M.CAPTURE.steps()),3); run_common(c); c.report()
