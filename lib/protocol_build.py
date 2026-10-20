"""Thin entry point used by data-backed named protocol pages."""
from __future__ import annotations

from pathlib import Path

from batch_page import render


def build(protocol, output: Path) -> None:
    output.write_text(render(protocol), encoding="utf-8")
    print(f"wrote {output}  ({output.stat().st_size:,} bytes)")
