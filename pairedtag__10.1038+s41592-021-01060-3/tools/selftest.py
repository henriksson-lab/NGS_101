#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import pairedtag as M
from checks import Check
c=Check(); c.section("PairedTag Supplementary Table 1 transcription")
c("pMENTs", M.PMENTS, "CTGTCTCTTATACACATCT")
c("AdaptorA", M.ADAPTOR_A, "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG")
c("Linker-R02", M.LINKER_R02, "CGAATGCTCTGGCCTCTCAAGCACGTGGAT")
c("Blocker-R02", M.BLOCKER_R02, "ATCCACGTGCTTGAGAGGCCAGAGCATTCG")
c("Linker-R03", M.LINKER_R03, "GGTCTGAGTTCGCACCGAAACATCGGCCAC")
c("Quencher-R03", M.QUENCHER_R03, "GTGGCCGATGTTTCGGTGCGAACTCAGACC")
c("Anchor-FokI-GH", M.ANCHOR_WRITTEN,
  "AAGCAGTGGTATCAACGCAGAGTGAAGGATGT"+"G"*9+"H")
c("P5c-NNDC-FokI", M.P5C_NNDC,
  "NNDCAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTG")
c("P7XX", M.P7XX,
  "CAAGCAGAAGACGGCATACGAGAT"+"N"*6+"GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC")
c("barcode lengths", (len(M.BC1_COMMON), M.R3_UMI_LENGTH,
                       M.R3_BARCODE_LENGTH, M.ROUND2.barcode_length),
  (15, 10, 7, 7))
c("published Read 2 barcode windows", M.PUBLISHED_READ2_BARCODE_WINDOWS,
  ((10, 13), (47, 50), (84, 87)))
c.report()
