#!/usr/bin/env python3
"""Build the FLASH-seq protocol schematic."""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHARED = HERE.parents[1] / "smart-seq__10.1038+nbt.2282" / "tools"
sys.path.insert(0, str(SHARED))

from page_parts import render_page  # noqa: E402

OUT = HERE.parent / "flash-seq.html"


def main() -> None:
    OUT.write_text(render_page("FLASH-seq"), encoding="utf-8")
    print(f"wrote {OUT}  ({OUT.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
