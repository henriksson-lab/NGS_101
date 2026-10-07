#!/usr/bin/env python3
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import splitseq  # noqa:F401,E402
from checks import Check  # noqa:E402
check=Check(); check.report()
