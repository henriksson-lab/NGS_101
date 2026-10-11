from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parents[1] / "lib")]
import chem
from batch_page import render

OUT = HERE.parent / "hicar.html"
if __name__ == "__main__":
    OUT.write_text(render(chem), encoding="utf-8")
    print(f"wrote {OUT}")
