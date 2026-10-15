#!/usr/bin/env python3
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from checks import Check
from chem import (BARCODE_1_AS_WRITTEN, GATE, GATC, RT_5N_3G, final_library)

c = Check()
c.section("Light-Seq source transcription and asymmetric library")
c("RT.5N.3G", RT_5N_3G, "TTTACACGATTGAGTTATNNNNNGGG")
c("GATE", GATE, "GGAGTTGGAGTGAGTGGATGAGTGATG")
c("GATC", GATC, "GAGAATGTGAGTGAAGATGTATGGTGA")
c("barcode 1 retains cnvK and dye notation", BARCODE_1_AS_WRITTEN,
  "GGAGTTGGAGTGAGTGGATGAGTGATGDDDDDDDDDDDDTATGGATGAGTTATATAACTCA[cnvK]TCGTGTAAAT[Cy5-3]")
lib, primers = final_library()
c("all four run primer roles land", tuple(p.role for p in primers),
  ("Read 1", "Index 1 (i7)", "Index 2 (i5)", "Read 2"))
c.report()
