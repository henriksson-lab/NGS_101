#!/usr/bin/env python3
"""Check source-printed Seq-Well S3 oligos against the original supplement."""
from __future__ import annotations

import sys
import os
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import seqwell_s3 as S
from checks import Check, Source, have

check = Check()
source_root = Path(os.environ.get("CHEM_DATA", HERE.parents[1] / "_data")) / "sources" / HERE.parent.name
original = Source(source_root / "j.immuni.2020.09.015_supp_mmc1.pdf.txt",
                  "original Document S1; reacquire with tools/get_sources.py")
if have(check, original, label="original oligo table"):
    text = original.path.read_text()
    # Read the original Sequences table rather than a copied model constant.
    table = text[text.rfind("Sequences"):]
    compact = re.sub(r"\s+", "", table)
    check("P5-SMART complete sequence and terminal phosphorothioates",
          S.P5_HYBRID[:-2] + "*A*C" in compact)
    check("source-printed modified TSO", "".join(s.top for s in S.tso_oligo()) in compact)
    check("source-printed S3 randomer", "".join(s.top for s in S.randomer()) in compact)
    check("source-printed custom Read 1", S.CUSTOM_R1 in compact)
check.report()
