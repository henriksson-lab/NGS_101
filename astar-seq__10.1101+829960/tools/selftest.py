#!/usr/bin/env python3
"""Checks for ASTAR source transcriptions not enforceable by the constructors.

Complementarity, strand placement and sequencing-primer landing are enforced while
building the page by Construct, Scene and seqprimers.section respectively.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import astar as A
from checks import Check
import nextera as nx

check = Check()

check.section("verbatim oligos from Supplementary Table 5")
check("C1-P2-T31", A.C1_DT,
      "GGCGACAACACCGATTGATCACGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")
check("C1-P2-RNA-TSO (bases; all are ribonucleotides)", A.C1_TSO,
      "GGCGACAACACCGATTGATCAGGG")
check("C1-P2-PCR-2", A.C1_PCR, "GGCGACAACACCGATTGATCA")
check("the two ATAC qPCR oligos",
      (A.ATAC_QPCR_S5, A.ATAC_QPCR_S7), (nx.ADAPTOR_S5, nx.ADAPTOR_S7))

check.report()
