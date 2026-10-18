#!/usr/bin/env python3
"""Construct native graph and page through molecular publication contracts."""
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parents[1]/'lib')]
import chem
from batch_page import render
render(chem)
print('PASS  Seq-Scope native graph and public-page construction')
