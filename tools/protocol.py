#!/usr/bin/env python3
"""Build or self-test one wrapper-free protocol module.

    python3 tools/protocol.py build <protocol-dir>
    python3 tools/protocol.py test <protocol-dir>
"""
from __future__ import annotations

import argparse
from pathlib import Path
import runpy
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))
from protocol_pages import load_page_module, page_module  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("build", "test"))
    parser.add_argument("directory", help="protocol directory name")
    args = parser.parse_args()
    directory = (ROOT / args.directory).resolve()
    if directory.parent != ROOT or not directory.is_dir():
        parser.error("directory must name one protocol under the repository root")
    declared = page_module(directory)
    if declared is None:
        parser.error("protocol has no PAGE_FILE module")
    module = load_page_module(directory)
    if args.action == "test":
        if callable(getattr(module, "selftest", None)):
            module.selftest()
        else:
            runpy.run_path(str(directory / "tools" / "selftest.py"), run_name="__main__")
    else:
        rendered = module.render_page()
        if not isinstance(rendered, str) or not rendered.strip():
            raise ValueError("render_page() must return nonempty HTML")
        out = directory / declared[1]
        out.write_text(rendered, encoding="utf-8")
        print(f"wrote {out}  ({out.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()
