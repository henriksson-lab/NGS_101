"""Small source/model boundary checks for one STRT protocol."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))

from checks import Check
import strt as S


def run(protocol: str) -> None:
    check = Check()
    check.section("published oligo assembled from model components")
    if protocol == S.STRT:
        check("STRT-V3-T30", S.STRT_DT, S.SOURCE_STRT_V3)
    elif protocol == S.C1:
        assembled = "".join(x.top for x in S.tn5_indexed_top("CGTCTAAT"))
        check("C1-TN5-1", assembled, S.SOURCE_C1_TN5_1)
    elif protocol == S.TWO_I:
        printed = S.TWO_I_WELL_PRIMER.replace("N" * S.WELL_INDEX_LEN, "X" * S.WELL_INDEX_LEN)
        check("DI-P1A-idx-P1B", printed, S.SOURCE_TWO_I_WELL)
    else:
        raise ValueError(protocol)
    check.report()
