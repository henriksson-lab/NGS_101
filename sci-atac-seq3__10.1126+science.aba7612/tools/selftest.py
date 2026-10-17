#!/usr/bin/env python3
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "sci-atac-seq__10.1126+science.aab1601" / "tools"))
sys.path.insert(0, str(HERE.parents[1] / "lib"))
sys.path.insert(0, str(HERE))
from checks import Check
from selftest_common import run
import sciatac
import sci3_graph

# Table S7, sheet 3lv2_splintandseqprimers (both splints are printed with /3ddC/).
SOURCE_N5_SPLINT = "GCCGACGACTGATTAC"
SOURCE_N7_SPLINT = "CACGAGACGACAAGTC"
SOURCE_IDX2 = "ACACATCTGACGCTGCCGACGACTGATTAC"
SOURCE_SEQ_R1 = "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
SOURCE_SEQ_R2 = "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"

if __name__ == "__main__":
    check = Check()
    check.section("Table S7 splints and Index 2 primer derived from the model")
    check("3LV2_N5_Splint", sci3_graph.splint_sequence("N5"), SOURCE_N5_SPLINT)
    check("3LV2_N7_Splint", sci3_graph.splint_sequence("N7"), SOURCE_N7_SPLINT)
    check("3lv2_seq_IDX2 reconstruction", sciatac.SCI3_I2, SOURCE_IDX2)
    primers = {p.role: p.seq for p in sciatac.SCI3_PRIMERS}
    check("3lv2_seq_R1 is the declared Read 1 primer", primers["Read 1"], SOURCE_SEQ_R1)
    check("3lv2_seq_R2 is the declared Read 2 primer", primers["Read 2"], SOURCE_SEQ_R2)
    check.report()
    run(sciatac.SCI3)
