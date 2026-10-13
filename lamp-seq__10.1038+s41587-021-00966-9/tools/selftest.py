#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
from checks import Check
from isothermal_protocols import (LAMP_FIP, LAMP_FIP_BARCODE,
                                  LAMP_POOL_FWD, LAMP_POOL_REV, get)
c=Check()
c("representative FIP barcode is embedded once", LAMP_FIP.count(LAMP_FIP_BARCODE), 1)
c("published pool-primer lengths", (len(LAMP_POOL_FWD), len(LAMP_POOL_REV)), (72, 68))
c("finished graph is publishable", get("lamp-seq").protocol_graph().publishable, True)
c.report()
