#!/usr/bin/env python3
# Construction and primer landing are enforced while chem.py builds its libraries.
import sys
from pathlib import Path
H = Path(__file__).resolve().parent
sys.path[:0] = [str(H), str(H.parents[1] / "lib")]
import chem
import storage
from checks import Check, Source, have

chem.protocol_graph()
check = Check()
check.section("Supplementary Data 1 Primers sheet transcription")
data1 = Source(storage.sources_root() / "redeem__10.1038+s41586-024-07066-z"
               / "s41586-024-07066-z_supp_41586_2024_7066_MOESM5_ESM.xlsx.txt",
               "Weng et al. 2024 Supplementary Data 1 "
               "(python3 tools/get_sources.py --slugs redeem__10.1038+s41586-024-07066-z)")
if have(check, data1):
    rows = dict(line.split("\t")[:2] for line in data1.path.read_text().splitlines()
                if line.startswith(("Nextera_N70", "P5\t", "P7\t")))
    check("N701-N706 primers equal P7 + index + S7 as printed",
          {k: chem.n70x_primer(v) for k, v in chem.N70X_INDEXES.items()},
          {k.removeprefix("Nextera_"): v for k, v in rows.items() if k.startswith("Nextera_")})
    check("post-capture P5/P7 as printed", (chem.REDEEM_P5, chem.REDEEM_P7),
          (rows["P5"], rows["P7"]))
check.report()
