#!/usr/bin/env python3
"""SMART-seq2 source transcription and inherited-linker evidence boundary."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]

import smartseq2 as s2
from checks import Check

check = Check()
# Hagemann-Jensen et al. 2020 Supplementary Table 1, condition 1: the
# Nature Protocols form used by later Smart-seq2 descendants.
check("Nature Protocols TSO, with ribo and LNA modifications", s2.TSO,
      "AAGCAGTGGTATCAACGCAGAGTACATrGrG+G")
check("Nature Protocols oligo-dT", s2.OLIGO_DT,
      "AAGCAGTGGTATCAACGCAGAGTAC" + "T" * 30 + "VN")

check("an inherited TSO linker is inferred by construction",
      s2.preamp_template().get("TSO linker").inferred)
check("a source-printed ACAT linker can be declared established",
      s2.preamp_template(linker_inferred=False).get("TSO linker").inferred,
      False)
check.report()
