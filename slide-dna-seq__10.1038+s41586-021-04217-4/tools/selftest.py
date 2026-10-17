#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import slide_dna_seq as P
from checks import Check
c=Check();c.section("slide-DNA source transcription");c("custom read primers",(P.READ1,P.READ2),("GCTTTGCTAACGGTCGAGAGATGTGTATAAGAGACAG","CGGATGTTGCACCAGCAGATGTGTATAAGAGACAG"))
c("custom index primers",(P.INDEX1_SEQ,P.INDEX2_SEQ),("GCTGGTGCAACATCCGTAAGCCTGGTGCCAAC","CTCGACCGTTAGCAAAGCTCCGCCCATGTGTC"))
c("bead oligo as printed",P.BEAD_OLIGO_5+"NNNNNNN"+P.BEAD_LINKER+"NNNNNNN"+P.BEAD_OLIGO_3,
  "GCCGCTACACGACGCTCTTCCGATCTNNNNNNNTCTTCAGCGTTCCCGAGANNNNNNNGCTCGGACACATGGGCG")
c.report()
