#!/usr/bin/env python3
"""Source transcription for the separately labelled 2019 comparator inventory."""
from pathlib import Path
import re
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import chem
from checks import Check, Source, have
from storage import sources_root

c = Check()
source = Source(sources_root() / HERE.parent.name / "PMC6615061_PMC6615061.xml.txt",
                "2019 comparator, reagent section: PCR primers 1.0 P5 and 2.0 P7")
if have(c, source, label="2019 comparator primer inventory"):
    text = source.path.read_text()
    for name, sequence in (("1.0 P5", chem.PCR_PRIMER_P5), ("2.0 P7", chem.PCR_PRIMER_P7)):
        match = re.search(r"PCR primer " + re.escape(name) + r"\s*\(5[’′']\s*([ACGT\s]+)3[’′']", text)
        c("2019 comparator " + name + " source transcription",
          sequence, re.sub(r"\s+", "", match.group(1)) if match else "")
c.report()
