from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from checks import Check
import chem
c = Check()
c("two known-length WGA4 adapters", [s.length_bp for s in chem.FINAL_LIBRARY if "WGA4 adapter" in s.name], [30, 30])
c("original run omits Read 2", [p.role for p in chem.SEQ_PRIMERS], ["Read 1", "Index 1 (i7)", "Index 2 (i5)"])
c.report()
