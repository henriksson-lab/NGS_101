#!/usr/bin/env python3
"""Primary-source transcription checks for 5-prime direct-capture Perturb-seq."""
from __future__ import annotations
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import direct_capture as C
import seqprimers as sp
from checks import Check, run_common

check = Check()
check.section("5-prime direct-capture Perturb-seq source transcription")
check("oJR160", C.OJR160,
      "AAGCAGTGGTATCAACGCAGAGTACCAAGTTGATAACGGACTAGCC")
check("oJR163", C.OJR163,
      "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT")
check("oJR165", C.OJR165,
      "CAAGCAGAAGACGGCATACGAGATAGGAGTCCGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGAGTACCAAGTTGATAACGGACTAGCC")
check("printed approximately 250-bp product", len(C.final_library()), 251)
check("all run primers land on the final library",
      sp.verify(C.final_library(), C.sequencing_primers(),
                required_roles=("Read 1", "Index 1 (i7)", "Read 2")), [])
run_common(check)
check.report()
