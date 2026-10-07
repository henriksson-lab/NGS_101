"""Source/model boundary checks for one SCRB protocol."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from checks import Check
import scrb as S


def run(protocol: str) -> None:
    check = Check()
    check.section("published oligo assembled from canonical components")
    if protocol == S.SCRB:
        check("E3V6NEXT", S.assembled_e3(), S.SOURCE_E3)
    elif protocol == S.MCSCRB:
        assembled = S.PCR_HANDLE + "GGG"
        check("unblocked E5V6NEXT", assembled, S.SOURCE_MC_TSO)
    else:
        raise ValueError(protocol)
    check.report()

