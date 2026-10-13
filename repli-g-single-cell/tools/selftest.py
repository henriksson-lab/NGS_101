#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
from checks import Check
from isothermal_protocols import get
c=Check(); c("finished graph is publishable", get("repli-g-single-cell").protocol_graph().publishable, True); c.report()
