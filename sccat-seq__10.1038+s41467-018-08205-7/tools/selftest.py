#!/usr/bin/env python3
"""Source-transcription checks for scCAT-seq."""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chem
import smartseq2 as ss2
from checks import Check
c=Check()
# Liu et al. 2019 Methods print the RNA-branch oligos verbatim, so the ACAT linker the
# shared SMART-seq2 model supplies is this paper's own published form, not an inheritance
# from either Picelli paper.  SMART-seq2 publishes two linkers (AC for the chosen rGrG+G
# oligo in Picelli 2013 Supplementary Table 2, ACAT in the Nature Protocols form); these
# checks pin the shared constants to what scCAT-seq itself printed.
c("published template-switching oligo",ss2.TSO,"AAGCAGTGGTATCAACGCAGAGTACATrGrG+G")
c("published cDNA PCR primer",ss2.ISPCR,"AAGCAGTGGTATCAACGCAGAGT")
c("published modified oligo-dT primer",ss2.OLIGO_DT,
  "AAGCAGTGGTATCAACGCAGAGTAC"+"T"*30+"VN")
c.report()
