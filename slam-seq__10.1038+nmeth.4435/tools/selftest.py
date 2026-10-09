#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent; sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import slam_seq as S
from checks import Check
c=Check(); c.section("source transcription")
c("standard alkylation conditions", S.ALKYLATION, {"IAA_mM":10,"DMSO_percent":50,"phosphate_mM":50,"pH":8,"temperature_C":50,"minutes":15})
c.report()
