from pathlib import Path
import sys
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
from checks import Check
import chem
c = Check()
c("GEM barcode length", len(next(s for s in chem.FINAL_LIBRARY if s.name.startswith("16-nt"))), 16)
c("random-primer UMI length", len(next(s for s in chem.FINAL_LIBRARY if s.name.startswith("8-nt"))), 8)
c.report()
