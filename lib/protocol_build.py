"""Thin entry point used by data-backed named protocol pages."""
from __future__ import annotations

from pathlib import Path

from batch_page import render
from protocol_specs import get


def build(protocol, output: Path) -> None:
    spec = get(protocol) if isinstance(protocol, str) else protocol
    output.write_text(render(spec), encoding="utf-8")
    print(f"wrote {output}  ({output.stat().st_size:,} bytes)")
