#!/usr/bin/env python3
"""Source transcription checks against the original selected NEB oligo manual."""
from __future__ import annotations
import os
import re
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import nebnext_user as N
from checks import Check, Source, have
check = Check()
bundle = Path(os.environ.get("CHEM_DATA", HERE.parents[1] / "_data")) / "sources" / HERE.parent.name
source = Source(bundle / "manuale7600.pdf.txt", "NEB E7600 v6.0_6/24; ref/MANIFEST.md")
if have(check, source, label="selected original NEB oligo manual"):
    original = re.sub(r"\s+", "", source.path.read_text())
    check("full modified hairpin in original E7600 p3",
          "/5Phos/" + N.hairpin_oligo().top().replace("U", "dU")[:-1] + "-s-T" in original)
    for primer in (N.i5_primer(), N.i7_primer()):
        sequence = primer.top()
        check(primer.name + " original complete oligo including terminal protection",
              sequence[:-1] + "*" + sequence[-1] in original)
second = Source(bundle / "manuale7780.pdf.txt", "NEB E7780 Set 2 original p3; independent same-adaptor source")
if have(check, second, label="second original hairpin source"):
    original = re.sub(r"\s+", "", second.path.read_text())
    check("full modified hairpin also printed by E7780 p3",
          "/5Phos/" + N.hairpin_oligo().top().replace("U", "dU")[:-1] + "-s-T" in original)
check.report()
