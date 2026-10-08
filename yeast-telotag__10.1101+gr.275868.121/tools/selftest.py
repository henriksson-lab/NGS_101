#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import telotag as M
from checks import Check,run_common
c=Check(); c.section("TeloTag source choices"); c("poly(A) represented as A only",set(M.POLYA),{"A"}); c("oligo(dT) pairs tail",len(M.POLYA),18); run_common(c); c.report()
