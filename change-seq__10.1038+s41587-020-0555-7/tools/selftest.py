#!/usr/bin/env python3
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from checks import Check
from chem import CHANGE_TRANSPOSON_BOTTOM, CHANGE_TRANSPOSON_TOP

c = Check()
c.section("CHANGE-seq Supplementary Protocol transcription")
c("oCRL225 molecular bases", CHANGE_TRANSPOSON_TOP, "ACGUAGATGTGTATAAGAGACAG")
c("oCRL226 molecular bases", CHANGE_TRANSPOSON_BOTTOM, "CTGTCTCTTATACACATCTACGT")
c.report()
