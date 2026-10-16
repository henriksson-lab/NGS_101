"""Locations of downloaded documents, text conversions and other untracked data."""
from __future__ import annotations

import os
from pathlib import Path

DEFAULT_DATA_ROOT = Path("/husky/henriksson/NGS101_papers")


def data_root() -> Path:
    """Use external storage by default; CHEM_DATA explicitly selects another root."""
    return Path(os.environ.get("CHEM_DATA", DEFAULT_DATA_ROOT)).expanduser()


def sources_root() -> Path:
    return data_root() / "sources"


def legacy_reference(relative: str | Path) -> Path:
    """An original formerly stored in a protocol ref/ directory or the pdf/ folder."""
    relative = Path(relative)
    if relative.is_absolute() or ".." in relative.parts:
        raise ValueError("legacy reference must be relative to the repository root")
    return data_root() / "legacy" / relative
