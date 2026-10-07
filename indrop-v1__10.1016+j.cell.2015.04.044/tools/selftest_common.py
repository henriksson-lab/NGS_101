"""Focused source/model checks for inDrop v1 and v2."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from checks import Check
import indrop as I


def run(protocol: str) -> None:
    check = Check(); check.section("published sequence assembled from model components")
    if protocol == I.V1:
        assembled = I.LEADER + I.T7_PROMOTER + I.MID + I.PE1[:24]
        check("shared acrydite primer bases (v2 supplement; v1 use remains inferred)", assembled, I.SOURCE_ACRYDITE_BASES)
    elif protocol == I.V2:
        check("PE2-N6", I.PE2 + "N" * 6, I.SOURCE_PE2_N6)
    else: raise ValueError(protocol)
    check.report()

