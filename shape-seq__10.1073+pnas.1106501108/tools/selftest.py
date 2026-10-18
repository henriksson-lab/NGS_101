#!/usr/bin/env python3
"""Build graph and public page to exercise molecular construction contracts."""
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import chem
from batch_page import render
render(chem)
print('PASS  SHAPE-Seq native graph and public-page construction')
