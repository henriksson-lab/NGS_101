#!/usr/bin/env python3
"""m6A-seq suite.

The final library is a role-only endpoint (no sequence segments or identifier
features), so no transcription check applies here; graph validity is enforced
when build_page.py renders the protocol graph.
"""
from pathlib import Path
import sys
H=Path(__file__).resolve().parent;sys.path[:0]=[str(H),str(H.parents[1]/"lib")]
from checks import Check
c=Check();c.report()
