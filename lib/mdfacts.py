"""
Computed facts inside Markdown notes, so prose cannot drift from the model.

    The probe is {{= crisprmip.PROBE_LEN }} nt.          -> the value, inline
    ```
    {{scene: crisprmip.probe_scene()}}                   -> a Scene's rows, as plain text
    {{oligo: crisprmip.probe_construct()}}               -> a Construct, 5'->3' with labels
    ```

Expressions are evaluated with every lib module and every protocol module (each
`<protocol>/tools/*.py` except build_page / selftest / show_*) importable by its own name,
plus `revcomp`, `tm`, `len`. A note inside a protocol directory resolves module names from
its own `<protocol>/tools/` first, so `chem` in a note means that protocol's `chem.py`.
A protocol module name defined in several protocol directories (such as `chem`) that is not
the note's own is ambiguous and an explicit error, never a silent pick of another protocol.
A failing expression is a build error naming the note and line -- a fact that no longer
resolves is exactly what this exists to catch.
"""

from __future__ import annotations

import importlib
import re
import sys
from pathlib import Path
from protocol_pages import page_module

ROOT = Path(__file__).resolve().parents[1]
FACT = re.compile(r"\{\{\s*(=|scene:|oligo:|duplex:)\s*(.+?)\s*\}\}")
_SKIP = {"build_page", "selftest"}
_base: dict | None = None
_dirs: dict[Path, dict] = {}          # tools dir -> {stem: module}
_stems: dict[str, list[Path]] = {}    # protocol module stem -> tools dirs defining it
_local: dict[Path, dict] = {}         # tools dir -> full namespace for its notes


class AmbiguousModule(Exception):
    pass


class _Ambiguous:
    """Stands for a protocol module name that several protocol directories define."""

    def __init__(self, name: str, dirs: list[Path]):
        self._name, self._dirs = name, dirs

    def _fail(self):
        where = ", ".join(d.parent.name for d in self._dirs[:3])
        more = f" and {len(self._dirs) - 3} more" if len(self._dirs) > 3 else ""
        raise AmbiguousModule(
            f"ambiguous module name {self._name!r}: defined in {len(self._dirs)} protocol "
            f"directories ({where}{more}); only a note inside one of them may use it")

    def __getattr__(self, attr):
        if attr.startswith("__"):
            raise AttributeError(attr)
        self._fail()

    def __call__(self, *args, **kwargs):
        self._fail()


def _modules(d: Path) -> list[str]:
    return [f.stem for f in sorted(d.glob("*.py"))
            if f.stem not in _SKIP and not f.stem.startswith(("show_", "_"))]


def _tool_dirs() -> list[Path]:
    return sorted(p for p in ROOT.glob("*/tools")
                  if (p / "build_page.py").exists() or (p / "catalogue.py").exists()
                  or page_module(p.parent) is not None)


def _load_dir(d: Path) -> dict:
    """Import one protocol tools dir; its own modules shadow same-named ones elsewhere."""
    if d in _dirs:
        return _dirs[d]
    names = _modules(d)
    shared = [n for n in names if len(_stems.get(n, ())) > 1 or n not in _stems]
    saved = {n: sys.modules.pop(n) for n in shared if n in sys.modules}
    sys.path.insert(0, str(d))
    found: dict = {}
    try:
        for name in names:
            try:
                found[name] = importlib.import_module(name)
            except Exception:          # a module that needs data we do not have
                pass
    finally:
        sys.path.remove(str(d))
        for n in shared:               # never leak this directory's `chem` to another
            sys.modules.pop(n, None)
        sys.modules.update(saved)
    _dirs[d] = found
    return found


def namespace() -> dict:
    """lib modules, helpers and every protocol module whose name is unambiguous."""
    global _base
    if _base is None:
        lib = ROOT / "lib"
        tools = _tool_dirs()
        for d in [lib, *tools]:
            if str(d) not in sys.path:
                sys.path.append(str(d))
        for d in tools:
            for name in _modules(d):
                _stems.setdefault(name, []).append(d)
        from chemdraw import revcomp, tm
        ns: dict = {"revcomp": revcomp, "tm": tm, "len": len, "round": round}
        for f in sorted(lib.glob("*.py")):
            name = f.stem
            if name in _SKIP or name.startswith(("show_", "_")) or name in ns:
                continue
            try:
                ns[name] = importlib.import_module(name)
            except Exception:          # a module that needs data we do not have
                pass
        for d in tools:
            for name, module in _load_dir(d).items():
                if name in ns:
                    continue
                if len(_stems[name]) > 1:
                    ns[name] = _Ambiguous(name, _stems[name])
                else:
                    ns[name] = module
        _base = ns
    return _base


def namespace_for(where: str = "") -> dict:
    """The namespace for a note at repository-relative path `where`."""
    base = namespace()
    parts = Path(where).parts if where else ()
    if len(parts) < 2:
        return base
    d = ROOT / parts[0] / "tools"
    if not d.is_dir():
        return base
    if d not in _local:
        ns = dict(base)
        lib = {f.stem for f in (ROOT / "lib").glob("*.py")}
        for name, module in _load_dir(d).items():
            if name not in lib:       # lib modules stay global
                ns[name] = module
        _local[d] = ns
    return _local[d]


def expand(md: str, where: str = "") -> str:
    def one(m: re.Match) -> str:
        kind, expr = m.group(1), m.group(2)
        line = md.count("\n", 0, m.start()) + 1
        try:
            val = eval(expr, namespace_for(where))      # noqa: S307 -- trusted repo code
        except Exception as e:
            raise ValueError(f"{where}:{line}: {{{{{kind} {expr}}}}} failed: {e!r}") from e
        if kind == "scene:":
            return "\n".join(r.plain().rstrip() for r in val.rows())
        if kind in ("oligo:", "duplex:"):
            from chemdraw import annotation_rows, strand_row
            rows = [strand_row(val, "top")]
            if kind == "duplex:":
                rows.append(strand_row(val, "bottom"))
            rows += annotation_rows(val)
            return "\n".join(r.plain().rstrip() for r in rows)
        return str(val)
    return FACT.sub(one, md)
