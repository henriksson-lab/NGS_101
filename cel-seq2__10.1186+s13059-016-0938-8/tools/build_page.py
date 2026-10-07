#!/usr/bin/env python3
"""Build the CEL-Seq2 protocol schematic."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHARED = HERE.parents[1] / "cel-seq__10.1016+j.celrep.2012.08.003" / "tools"
sys.path.insert(0, str(SHARED))

from page_parts import render_page

OUT = HERE.parent / "cel-seq2.html"


def main() -> None:
    OUT.write_text(render_page("CEL-Seq2"), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
