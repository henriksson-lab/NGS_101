#!/usr/bin/env python3
"""Check source-distinct assay stages that the role-level model must retain."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import superrca as M
from checks import Check

check = Check()
w = M.workflow()
check("each reaction consumes the preceding output",
      all(a.before is b.after for a, b in zip(w.reactions[1:], w.reactions[:-1])), True)
check("two distinct RCA reactions surround padlock genotyping",
      [r.action for r in w.reactions if "Phi29" in r.action],
      ["Phi29 copies the target circle from the ligation-template primer",
       "Added primers and Phi29 perform secondary RCA"])
check("amplicon circularization is splinted", M.CIRCULARIZATION.kind,
      "splinted single-strand ligation")
check.report()
