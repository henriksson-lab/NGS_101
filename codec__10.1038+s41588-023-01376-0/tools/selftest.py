#!/usr/bin/env python3
"""Compare the representative CODEC oligos with original Supplementary Table 5."""
import os
import re
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))
import chem as C
from checks import Check, Source, have
check = Check()
root = Path(os.environ.get("CHEM_DATA", HERE.parents[1] / "_data"))
src = Source(root / "sources" / HERE.parent.name / "s41588-023-01376-0_supp_41588_2023_1376_MOESM1_ESM.pdf.txt",
             "Original CODEC Supplementary Table 5; tools/get_sources.py --slugs " + HERE.parent.name)
if have(check, src, label="original Supplementary Table 5"):
    text = src.path.read_text()
    text = text[text.index("Supplementary Table 5."):]
    compact = re.sub(r"\s+", "", text)
    for name, seq in (("LD4-adap5-1", C.CODEC_ADAP5_1[:-1] + "*T"),
                      ("LD4-adap7-1", C.CODEC_ADAP7_1[:-1] + "*T"),
                      ("LD4-brid5-1", "P-" + C.CODEC_BRID5_1),
                      ("LD4-brid7-1", "P-" + C.CODEC_BRID7_1)):
        check(name + " original complete modified ordering sequence", name + seq in compact)
check.report()
