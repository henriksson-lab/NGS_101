#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from checks import Check
from chem import (CHEC_INDEXED_ADAPTER, CHEC_P5_PCR, CHEC_P7_PCR,
                                 CHEC_UNIVERSAL_ADAPTER)

c = Check()
c.section("cited low-input paired-end library oligos")
c("Universal Adapter bases", CHEC_UNIVERSAL_ADAPTER,
  "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT")
c("Indexed Adapter bases", CHEC_INDEXED_ADAPTER,
  "GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG")
c("P5 enrichment primer", CHEC_P5_PCR, "AATGATACGGCGACCACCGAG")
c("P7 enrichment primer", CHEC_P7_PCR, "CAAGCAGAAGACGGCATACGAG")
c.report()
