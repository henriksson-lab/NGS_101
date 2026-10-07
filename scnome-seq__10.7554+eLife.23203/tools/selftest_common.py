"""Source/model boundary checks for scNOMe-seq and scCOOL-seq."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from checks import Check
import nomecool as N


def run(protocol: str) -> None:
    check = Check(); check.section("published oligo assembled from canonical components")
    if protocol == N.SCNOME:
        check("forward_p5", N.FORWARD_P5, N.SOURCE_FORWARD_P5)
    elif protocol == N.SCCOOL:
        check("random primer 1", "".join(x.top for x in N.primer1_segments()), N.SOURCE_COOL_P1)
    else: raise ValueError(protocol)
    check.report()

