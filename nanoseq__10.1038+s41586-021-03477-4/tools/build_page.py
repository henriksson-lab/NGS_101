from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
from batch_page import render
import chem
OUT=HERE.parent/"nanoseq.html"
if __name__=="__main__":
 OUT.write_text(render(chem.SPEC),encoding="utf-8")
 print(f"wrote {OUT}")
