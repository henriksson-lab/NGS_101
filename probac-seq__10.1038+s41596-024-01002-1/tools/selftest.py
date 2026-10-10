#!/usr/bin/env python3
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
import probac_seq as P
from checks import Check
c=Check();c.section("ProBac identifier roles");c("probe feature is not cell identity",[(s.name,s.feature.role) for s in P.FINAL_LIBRARY if s.feature],[('cell barcode','cell_barcode'),('UMI','umi'),('gene-specific probe identity','feature_barcode'),('i7 reverse complement','sample_index')]);c.report()
