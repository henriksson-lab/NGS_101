from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
from ligation_free_conformation import get
globals().update(vars(get("gam")))
