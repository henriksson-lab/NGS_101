#!/usr/bin/env python3
"""CEL-Seq2 source transcription and construct boundary checks."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE),
                str(HERE.parents[1] / "cel-seq__10.1016+j.celrep.2012.08.003" / "tools"),
                str(HERE.parents[1] / "lib")]

import celseq as C
import chem
from checks import Check

check = Check()
# Hashimshony et al. 2016 Table S2 prints the CEL-Seq2 RT-primer design.
check("CEL-Seq2 RT-primer design", C.rt_primer(chem.PROTOCOL, chem.DESIGN).top(),
      "GCCGG" + C.T7_PROMOTER + "AGTTCTACAGTCCGACGATC" + "N" * 12
      + "T" * 24 + "V")
check("randomhexRT ordering string",
      chem.DESIGN.randomhex_tail + "N" * chem.DESIGN.randomhex_nt,
      "GCCTTGGCACCCGAGAATTCCANNNNNN")
geometry = chem.library_pcr_geometry()
check("vendor RP1/RPI1 face inward on CEL-Seq2 cDNA",
      (geometry.forward_hit.strand, geometry.reverse_hit.strand,
       geometry.forward_hit.end - geometry.forward_hit.start,
       geometry.reverse_hit.end - geometry.reverse_hit.start),
      ("bottom", "top", 16, 21))
check.report()
