#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
import telotag as M
from checks import Check,run_common
c=Check(); c.section("TeloTag source transcription"); c("published TeloTag core",M.TELOTAG_CORE,"CAGTCTACACATATTCTCTG"); c("five published primer variants",tuple(p[-3:] for p in M.TELOTAG_PRIMERS),("TTT","CAC","ACA","CCA","ACC")); c("printed oligo(dT) lengths",tuple(len(p)-len(M.TELOTAG_CORE)-(0 if i==0 else 3) for i,p in enumerate(M.TELOTAG_PRIMERS)),(35,36,36,36,36)); run_common(c); c.report()
