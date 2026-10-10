#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import direct_seq as P
from checks import Check
c=Check();c.section("Direct-seq source transcription");c("8A8G tract",P.CAPTURE,"AAAAAAAAGAAAAAAAGAAAAAAAGAAAAA");c.report()
