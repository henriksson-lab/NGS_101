#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "lib"))
from batch_ngs import seg, truseq_full_length_adapter, truseq_library
from checks import Check
import chem

c = Check()
# Construction validates the canonical library and all declared primer landing sites.
m = chem
c("PCR-free workflow has no PCR reaction", any("PCR" in r.action for r in m.workflow().reactions), False)

# Source transcription: which canonical oligo goes where in the ligated full-length
# adaptor. Illumina Adapter Sequences #1000000002694 v22 prints the two arms of the
# TruSeq dual-index adaptor and the index-less universal adaptor; the index windows are
# the only part the model keeps as a placeholder.
DUAL = truseq_full_length_adapter()
c("full-length i5 arm as printed (Index 2 (i5) Adapters)",
  DUAL.universal_sequence().replace("J" * 8, "[i5]"),
  "AATGATACGGCGACCACCGAGATCTACAC[i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT")
c("full-length i7 arm as printed (Index 1 (i7) Adapters)",
  DUAL.indexed_sequence().replace("I" * 8, "[i7]"),
  "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[i7]ATCTCGTATGCCGTCTTCTGCTTG")
c("single-index universal arm as printed (TruSeq Universal Adapter)",
  truseq_full_length_adapter(dual_index=False).universal_sequence(),
  "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT")

# The PCR-free constructor differs from the amplified one in provenance, not in bases:
# it must not invent or drop a base relative to the established TruSeq shell.
amplified, _ = truseq_library([seg("genomic DNA insert", "N" * 42, placeholder=True)],
                              "amplified comparison")
c("PCR-free library bases equal the canonical TruSeq shell",
  m.FINAL_LIBRARY.top(), amplified.top())
c.report()
