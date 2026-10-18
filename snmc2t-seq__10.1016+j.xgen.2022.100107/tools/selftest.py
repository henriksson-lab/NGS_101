#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import snmc2t_seq as P
from checks import Check
c=Check();c.section("Published2022 modified oligo transcription")
c("dT30VN_4 DNA deoxyuridines",P.DT30VN,"AAGCAGUGGUAUCAACGCAGAGUACUTTTTTUTTTTTUTTTTTUTTTTTUTTTTTVN")
c("N6_2",P.N6,"AAGCAGUGGUAUCAACGCAGAGUACNNNNNN")
c("TSO_3 RNA/LNA tail",P.TSO,"AAGCAGUGGUAUCAACGCAGAGUGAAUrGrG+G")
c("ISPCR23_2 deoxyuridines",P.ISPCR23,"AAGCAGUGGUAUCAACGCAGAGU")
c.report()
