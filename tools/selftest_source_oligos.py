#!/usr/bin/env python3
"""Local source transcription check for exact graph oligo inputs.

Run with local $CHEM_DATA/sources/<protocol-dir>/ text twins. Missing source
bundles SKIP. Protocol modules may declare SOURCE_OLIGO_EXCEPTIONS as a mapping
from graph material ID to the reason its sequence cannot occur verbatim.
"""
from __future__ import annotations

import argparse
import ast
import csv
import importlib.util
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))
sys.path.insert(0, str(ROOT / "tools"))
from storage import data_root
from source_oligo_check import check_graph, unrepresented_input_ids  # noqa: E402


def graph_module_paths(directory: Path) -> list[Path]:
    result = []
    for path in sorted((directory / "tools").glob("*.py")):
        if path.name == "selftest.py":
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        if any(isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
               and node.name == "protocol_graph" and not node.args.posonlyargs
               and not node.args.args and not node.args.kwonlyargs
               for node in tree.body):
            result.append(path)
    return result


def load(path: Path):
    sys.path.insert(0, str(path.parent))
    name = "source_check_" + path.parent.parent.name.replace("-", "_").replace("+", "_")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        sys.path.pop(0)
    return module


def clear_protocol_imports() -> None:
    """Prevent a local `chem` or `page_parts` import leaking into the next protocol."""
    for key, module in tuple(sys.modules.items()):
        filename = getattr(module, "__file__", None)
        if not filename:
            continue
        path = Path(filename).resolve()
        if path.parent.name == "tools" and path.parent.parent.parent == ROOT:
            sys.modules.pop(key, None)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("dirs", nargs="*", help="protocol directory names; default: every catalogue entry")
    parser.add_argument("--details", action="store_true",
                        help="show graph input IDs that lack an exact whole-strand state")
    args = parser.parse_args()
    with (ROOT / "catalogue" / "ours.tsv").open(encoding="utf-8", newline="") as handle:
        known = {row["dir"]: row for row in csv.DictReader(handle, delimiter="\t")}
    names = args.dirs or sorted(known)
    bad = [name for name in names if name not in known]
    if bad:
        parser.error("unknown protocol directory: " + ", ".join(bad))
    source_root = data_root() / "sources"
    failures = skips = checked = graphs = unrepresented = 0
    for name in names:
        for path in graph_module_paths(ROOT / name):
            try:
                module = load(path)
                sys.path.insert(0, str(path.parent))
                try:
                    graph = module.protocol_graph()
                finally:
                    sys.path.pop(0)
                missing = unrepresented_input_ids(graph)
                unrepresented += len(missing)
                if args.details and missing:
                    print(f"UNREPRESENTED {name}: {', '.join(missing)}")
                problems, omitted, count = check_graph(
                    graph, source_root / name,
                    getattr(module, "SOURCE_OLIGO_EXCEPTIONS", {}), known[name]["doi"])
            except Exception as exc:
                problems, omitted, count = [f"{name}/{path.name}: {type(exc).__name__}: {exc}"], [], 0
            finally:
                clear_protocol_imports()
            graphs += 1
            checked += count
            failures += len(problems)
            skips += len(omitted)
            for issue in problems:
                print("FAIL", issue)
            for issue in omitted:
                print("SKIP", issue)
    print(f"Source oligos: {checked} checked, {failures} failed, {skips} source bundles skipped, "
          f"{unrepresented} inputs without exact states, {graphs} graphs")
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
