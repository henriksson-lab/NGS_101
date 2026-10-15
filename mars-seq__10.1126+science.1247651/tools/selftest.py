#!/usr/bin/env python3
"""Source-transcription checks for MARS-seq (2014)."""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import marsseq as M  # noqa: E402
from checks import Check  # noqa: E402
check = Check()

# Jaitin et al. 2014, Science 343:776, Supplementary Material, Table S7 "Primers used
# during MARS-seq library construction" -- the publisher's own PDF, on disk as
# scg_jaitin-sm.pdf.  These three primers are assembled from illumina constants, so a
# change there would silently alter what this page presents as verbatim.  The scg_lib_structs
# page is not evidence for any of them and the oligo it calls MARS-seq2.0's RT2
# (CTACACGACGCTCTTCCGATCT) is deliberately absent from the model.
check("Table S7 barcoded RT primer (XXXXXX/NNNN modelled as N)", M.MARS1_RT1,
      "CGATTGAGGCCGGTAATACGACTCACTATAGGGGCGACGTGTGCTCTTCCGATCT"
      + "N" * 10 + "T" * 20 + "N")
check("Table S7 Second RT primer", M.MARS1_RT2_PRINTED, "TCTAGCCTTCTCGCAGCACATC")
check("Table S7 P5_Rd1 PCR forward", M.MARS1_P5_PCR,
      "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT")
check("Table S7 P7_Rd2 PCR reverse", M.MARS1_P7_PCR,
      "CAAGCAGAAGACGGCATACGAGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT")
# Table S8, row lig_NNNX4_ix1: /5Phos/GACTNNNAGATCGGAAGAGCGTCGTGTAG/3SpC3/
check("Table S8 ligation-adapter constant", M.LIG_CONSTANT, "AGATCGGAAGAGCGTCGTGTAG")
check.report()
