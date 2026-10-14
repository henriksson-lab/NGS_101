#!/usr/bin/env python3
"""Source-boundary guard for the LIANTI model.

The paper prints no oligo and its Supplementary Materials and Methods is unavailable, so
no LIANTI-specific base may be asserted.  The only check worth keeping is the one the
construction functions cannot make themselves: that the transposon's unpublished
spacer never acquires concrete bases again.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import lianti as L
from checks import Check

check = Check()
check.section("no unpublished LIANTI base is asserted")
spacers = [s for s in L.transposon_segments() if "spacer" in s.name]
check("spacer segment is a placeholder", [(s.placeholder, s.inferred) for s in spacers],
      [(True, True)])
check("spacer carries no concrete base", set("".join(s.top for s in spacers)), {"N"})
check.report()
