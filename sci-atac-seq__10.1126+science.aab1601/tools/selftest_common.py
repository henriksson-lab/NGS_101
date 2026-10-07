"""Source/model boundary checks for one sci-ATAC protocol."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from checks import Check
import illumina as il
import nextera as nx
import sciatac as S


def run(protocol: str) -> None:
    check = Check()
    check.section("published oligo assembled from canonical components")
    if protocol == S.SCI15:
        check("T5_1", "".join(x.top for x in S.t5("TATAGCCT")), S.SOURCE_T5_1)
    elif protocol == S.SCI18:
        assembled = il.P5 + "CTCCATCGAG" + nx.S5
        check("P5_1_PCR_Primer", assembled, S.SOURCE_2018_P5_1)
    elif protocol == S.SCI3:
        check("3LV2 Index 1 reconstruction", S.SCI3_I1, S.SOURCE_SCI3_I1)
    else:
        raise ValueError(protocol)
    check.report()

