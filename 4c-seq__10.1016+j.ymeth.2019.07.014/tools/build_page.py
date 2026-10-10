import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "lib"))
import chem
from protocol_build import build

OUT = HERE.parent / "4c-seq.html"
build(chem.SPEC, OUT)
