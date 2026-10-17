#!/usr/bin/env python3
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent; sys.path[:0]=[str(HERE),str(HERE.parents[1]/"lib")]
from checks import Check
from lamp_seq import LAMP_FIP, LAMP_FIP_BARCODE, LAMP_POOL_FWD, LAMP_POOL_REV
import lamp_seq as model
c=Check()
c("representative FIP barcode is embedded once", LAMP_FIP.count(LAMP_FIP_BARCODE), 1)
c("Supplementary Data 1 Pool_fwd_1/Pool_rev_1 lengths", (len(LAMP_POOL_FWD), len(LAMP_POOL_REV)), (62, 58))
c("finished graph is publishable", model.protocol_graph().publishable, True)
c.report()
