#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import chip_seq
from checks import Check, Source, have
from storage import data_root
c = Check()
# No oligo transcription check: no inspected original prints PCR primer 2.1 or the
# Adapter Oligo Mix bases (round-6 audit); the former constant-restating check was
# removed rather than presented as source verification. Construction is the check.
chip_seq.protocol_graph()
# Source transcription: the read length is stated only in the authors' Supplementary
# Data 1 README (Springer MOESM306), lines 18-21.
readme = Source(data_root() / "sources" / H.parent.name / "nmeth1068_41592_2007_BFnmeth1068_MOESM306_ESM_README.txt",
                "Robertson et al. 2007 Supplementary Data 1 zip; reacquire via tools/get_sources.py --slugs")
if have(c, readme, label="Supplementary Data 1 README read length"):
    text = " ".join(readme.path.read_text(errors="replace").split())
    c("README states the modelled read length",
      f"number of {chip_seq.FINAL_READ_LENGTH} bp Illumina reads" in text)
c.report()
