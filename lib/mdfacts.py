"""
Computed facts inside Markdown notes, so prose cannot drift from the model.

    The probe is {{= crisprmip.PROBE_LEN }} nt.          -> the value, inline
    ```
    {{scene: crisprmip.probe_scene()}}                   -> a Scene's rows, as plain text
    {{oligo: crisprmip.probe_construct()}}               -> a Construct, 5'->3' with labels
    ```

Expressions are evaluated with every lib module and every protocol module (each
`<protocol>/tools/*.py` except build_page / selftest / show_*) importable by its own name,
plus `revcomp`, `tm`, `len`. A failing expression is a build error naming the note and
line -- a fact that no longer resolves is exactly what this exists to catch.
"""

from __future__ import annotations

import importlib
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FACT = re.compile(r"\{\{\s*(=|scene:|oligo:|duplex:)\s*(.+?)\s*\}\}")
_SKIP = {"build_page", "selftest"}
_ns: dict | None = None


def namespace() -> dict:
    global _ns
    if _ns is None:
        lib = ROOT / "lib"
        tools = sorted(p for p in ROOT.glob("*/tools")
                       if (p / "build_page.py").exists() or (p / "catalogue.py").exists())
        for d in [lib, *tools]:
            if str(d) not in sys.path:
                sys.path.append(str(d))
        from chemdraw import revcomp, tm
        ns: dict = {"revcomp": revcomp, "tm": tm, "len": len, "round": round}
        for d in [lib, *tools]:
            for f in sorted(d.glob("*.py")):
                name = f.stem
                if name in _SKIP or name.startswith(("show_", "_")) or name in ns:
                    continue
                try:
                    ns[name] = importlib.import_module(name)
                except Exception:          # a module that needs data we do not have
                    pass
        _ns = ns
    return _ns


def expand(md: str, where: str = "") -> str:
    def one(m: re.Match) -> str:
        kind, expr = m.group(1), m.group(2)
        line = md.count("\n", 0, m.start()) + 1
        try:
            val = eval(expr, namespace())      # noqa: S307 -- trusted repo code
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
