#!/usr/bin/env python3
"""Source-transcription check for the CG000109 gel-bead TSO."""
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import tenx_5p_v1 as V
from checks import Check
check=Check()
got="".join(s.top for s in V.bead_tso()).replace("B","N").replace("U","N")
want="CTACACGACGCTCTTCCGATCT"+"N"*16+"N"*10+"TTTCTTATATGGG"
check("CG000109 gel-bead TSO transcription", got, want)
# CG000109 Rev D Table 2B (5' GEX), Protocol 4.1 / 6.3 / 6.5 as printed.
check("CG000109 Table 2B cDNA forward primer", V.CDNA_FORWARD, "CTACACGACGCTCTTCCGATCT")
check("CG000109 Table 2B cDNA reverse primer", V.CDNA_REVERSE, "AAGCAGTGGTATCAACGCAG")
check("CG000109 Table 2A cDNA reverse primer", V.CDNA_REVERSE_TABLE_2A, "AAGCAGTGGTATCAACGCAGAG")
check("CG000109 Table 2B RT primer", "".join(s.top for s in V.rt_primer()), "AAGCAGTGGTATCAACGCAGAGTAC"+"T"*30+"VN")
check("CG000109 Table 2B Read 2 adaptor", ("".join(s.top for s in V.adapter_top()), "".join(s.top for s in V.adapter_bottom())[::-1]), ("GATCGGAAGAGCACACGTCTGAACTCCAGTCAC", "TCTAGCCTTCTCG"))
check("CG000109 Table 2B P5 index primer", V.INDEX_FORWARD, "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC")
check("CG000109 Table 2B P7 index primer", V.INDEX_REVERSE, "CAAGCAGAAGACGGCATACGAGAT"+"N"*8+"GTGACTGGAGTTCAGACGTGT")
check.report()
