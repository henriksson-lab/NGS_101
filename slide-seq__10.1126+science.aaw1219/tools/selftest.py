#!/usr/bin/env python3
"""Published Table S1 decoding-primer transcription; construction is checked by the build."""
import re
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from storage import data_root
from checks import Check, Source, have
import chem

check = Check()
supp = Source(data_root() / "sources" / HERE.parent.name /
              "science.aaw1219_supp_NIHMS1036108-supplement-1.docx.txt",
              "PMC6927209 NIHMS1036108 supplement-1.docx (published author-manuscript supplement)")
if have(check, supp, label="Published Table S1"):
    rows = {re.sub(r"\s+", " ", m.group(1)).strip(): m.group(2)
            for m in re.finditer(r"^(.+?)\s*\t/5Phos/([ACGT]+)", supp.path.read_text(), re.M)}
    check("ten decoding primers transcribed as printed in Table S1",
          {name: rows.get(name) for name, _ in chem.DECODING_PRIMERS}, dict(chem.DECODING_PRIMERS))
check.report()
