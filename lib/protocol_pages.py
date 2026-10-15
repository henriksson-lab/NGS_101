"""Find and load protocol modules that build their own page without local wrappers."""
from __future__ import annotations

import ast
import importlib.util
from pathlib import Path
import re
import sys


def _declarations(path: Path) -> tuple[str | None, set[str]]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    page = None
    functions = set()
    for node in tree.body:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            functions.add(node.name)
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and (
                isinstance(node.targets[0], ast.Name) and
                node.targets[0].id == "PAGE_FILE" and
                isinstance(node.value, ast.Constant) and
                isinstance(node.value.value, str)):
            page = node.value.value
    return page, functions


def page_module(directory: Path) -> tuple[Path, str] | None:
    """Return the one local module with a literal PAGE_FILE declaration, if any."""
    found = []
    for path in sorted((directory / "tools").glob("*.py")):
        if path.stem in {"build_page", "selftest"}:
            continue
        page, functions = _declarations(path)
        if page is None:
            continue
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*\.html", page):
            raise ValueError(f"{path}: PAGE_FILE must be a plain HTML filename")
        if "render_page" not in functions:
            raise ValueError(f"{path}: page module needs render_page()")
        if "selftest" not in functions and not (directory / "tools" / "selftest.py").exists():
            raise ValueError(f"{path}: page module needs selftest() or local selftest.py")
        found.append((path, page))
    if len(found) > 1:
        raise ValueError(f"{directory}: multiple PAGE_FILE declarations")
    return found[0] if found else None


def load_page_module(directory: Path):
    declared = page_module(directory)
    if declared is None:
        raise ValueError(f"{directory}: no page module")
    path, _ = declared
    name = path.stem  # one protocol per subprocess; inherited modules see this instance
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ValueError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    # Source-owned models may import inherited molecular constructors by name.
    sys.path.insert(0, str(path.parent))
    spec.loader.exec_module(module)
    return module
