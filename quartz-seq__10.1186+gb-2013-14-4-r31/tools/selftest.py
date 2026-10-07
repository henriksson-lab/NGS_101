#!/usr/bin/env python3
"""Source-transcription checks for Quartz-Seq oligos."""
from __future__ import annotations
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "lib"))
import quartz_seq as Q
from checks import Check
check = Check(); check.section("2013 Table 1 transcription")
check("RT primer", Q.RT_PRIMER, "TATAGAATTCGCGGCCGCTCGCGATAATACGACTCACTATAGGGCGTTTTTTTTTTTTTTTTTTTTTTTT")
check("Tagging primer", Q.TAGGING_PRIMER, "TATAGAATTCGCGGCCGCTCGCGATTTTTTTTTTTTTTTTTTTTTTTT")
check("Suppression primer", Q.SUPPRESSION_PRIMER, "GTATAGAATTCGCGGCCGCTCGCGAT")
check("TRSU", Q.TRSU, "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT")
check("TRSI-2", Q.TRSI2, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG")
check.report()
