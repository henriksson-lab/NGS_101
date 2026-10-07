#!/usr/bin/env python3
"""Build the CEL-Seq protocol schematic."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from page_parts import render_page

OUT = HERE.parent / "cel-seq.html"


def main() -> None:
    OUT.write_text(render_page("CEL-Seq"), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
