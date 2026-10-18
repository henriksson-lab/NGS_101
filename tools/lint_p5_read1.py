#!/usr/bin/env python3
"""
Lint built molecular models for a malformed P5 / TruSeq Read 1 junction.

``illumina.P5`` ends in ACAC and ``illumina.TRUSEQ_READ1`` begins with ACAC.  In a
dual-indexed primer the i5 window sits between them, so both copies are real.  Without
an i5 window the two are one primer sequence, ``illumina.TRUSEQ_P5_FULL``
(...GAGATCTACACTCTTTCCC...), and the ACAC occurs once.  Placing the two constants as
adjacent segments gives ...GAGATCTACACACACTCTTTCCC...; dropping both copies gives
...GAGATCTCTCTTTCCC....  ``seqprimers`` still finds Read 1 in either case, so nothing
else catches it.

This lint imports every protocol module (each directory in its own subprocess, run in
parallel), collects every ``Construct`` it can reach -- module-level objects and the
states and sequencing libraries of ``protocol_graph()`` -- and flags any whose sequence,
on either strand, joins the P5 core to the Read 1 arm with other than exactly one ACAC.
It reports the offending segment boundary.  It checks sequence arithmetic only; whether
an i5 window should be present is a source question for the protocol.

A protocol whose authoritative source itself prints a primer with the duplicated (or
missing) ACAC declares it in the module that builds the construct:

    P5_READ1_JUNCTION_EXCEPTIONS = {"<construct name>": "<reason with source locator>"}

Matching findings are reported as ``exempted`` (not hidden) and do not fail the lint.
An exemption whose reason is empty, or that matches no finding, is itself a finding.

Run:  python3 tools/lint_p5_read1.py [protocol-dir ...]   (default: every protocol)
Exit status 1 if anything is found; directories whose modules fail to import are
listed as skipped, not as findings.
"""
from __future__ import annotations

import argparse
import concurrent.futures as cf
import importlib.util
import inspect
import json
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LIB = ROOT / "lib"
SKIP_STEMS = {"build_page", "selftest"}


def _patterns():
    sys.path.insert(0, str(LIB))
    import illumina as il
    from chemdraw import revcomp
    core = il.P5[:-4]
    if core + il.TRUSEQ_READ1 != il.TRUSEQ_P5_FULL or not il.P5.endswith("ACAC"):
        raise AssertionError("P5 / Read 1 overlap assumption no longer holds")
    arm = il.TRUSEQ_READ1[4:]
    forward = re.compile(core + "((?:ACAC)*)" + arm)
    reverse = re.compile(revcomp(arm) + "((?:GTGT)*)" + revcomp(core))
    return forward, reverse


def _constructs(module):
    """Every Construct reachable from module globals and its protocol graph."""
    from chemdraw import Construct
    seen, out = set(), []

    def add(obj, where, depth=0):
        if id(obj) in seen or depth > 3:
            return
        seen.add(id(obj))
        if isinstance(obj, Construct):
            out.append((where, obj))
        elif isinstance(obj, (list, tuple, set, frozenset)):
            for i, x in enumerate(obj):
                add(x, f"{where}[{i}]", depth + 1)
        elif isinstance(obj, dict):
            for k, x in obj.items():
                add(x, f"{where}[{k!r}]", depth + 1)

    for name, obj in vars(module).items():
        if not name.startswith("__"):
            add(obj, name)
    factory = getattr(module, "protocol_graph", None)
    if callable(factory) and _no_required_args(factory):
        graph = factory()
        for mid, material in getattr(graph, "materials", {}).items():
            add(material.state, f"protocol_graph().materials[{mid!r}].state")
        for mid, plan in getattr(graph, "sequencing", {}).items():
            add(plan.library, f"protocol_graph().sequencing[{mid!r}].library")
    return out


def _no_required_args(fn) -> bool:
    try:
        params = inspect.signature(fn).parameters.values()
    except (TypeError, ValueError):
        return False
    return not any(p.default is p.empty and p.kind in (p.POSITIONAL_ONLY,
                   p.POSITIONAL_OR_KEYWORD, p.KEYWORD_ONLY) for p in params)


def _boundary(con, start: int, end: int) -> str:
    """Name the segments covering the matched junction region."""
    names, pos = [], 0
    for s in con:
        after = pos + len(s)
        if after > start and pos < end and s.name not in names:
            names.append(s.name or "unnamed")
        pos = after
    return " | ".join(names)


def scan_directory(directory: Path) -> dict:
    """Worker: import one protocol's modules and return findings as a dict."""
    forward, reverse = _patterns()
    tools = directory / "tools"
    sys.path[:0] = [str(tools), str(LIB), str(ROOT)]
    findings, skipped, exempted = [], [], []
    for path in sorted(tools.glob("*.py")):
        if path.stem in SKIP_STEMS:
            continue
        try:
            spec = importlib.util.spec_from_file_location(path.stem, path)
            module = importlib.util.module_from_spec(spec)
            sys.modules[path.stem] = module
            spec.loader.exec_module(module)
            constructs = _constructs(module)
        except BaseException as exc:  # noqa: BLE001 - a broken module is reported, not fatal
            skipped.append(f"{path.name}: {type(exc).__name__}: {str(exc)[:160]}")
            continue
        exceptions = getattr(module, "P5_READ1_JUNCTION_EXCEPTIONS", {}) or {}
        if not isinstance(exceptions, dict):
            findings.append({"module": path.name, "where": "P5_READ1_JUNCTION_EXCEPTIONS",
                             "construct": "", "strand": "-", "segments": "-",
                             "problem": "exception declaration is not a dict"})
            exceptions = {}
        used = set()
        reported = set()
        for where, con in constructs:
            seq = con.top()
            for strand, pattern, unit in (("top", forward, "ACAC"), ("bottom", reverse, "GTGT")):
                for m in pattern.finditer(seq):
                    copies = len(m.group(1)) // 4
                    if copies == 1:
                        continue
                    key = (con.name, seq, m.start(), strand)
                    if key in reported:
                        continue
                    reported.add(key)
                    kind = "duplicated ACAC" if copies > 1 else "missing ACAC"
                    finding = {
                        "module": path.name, "where": where, "construct": con.name,
                        "strand": strand, "problem": f"{kind} ({copies} copies)",
                        "segments": _boundary(con, m.start(), m.end()),
                    }
                    reason = exceptions.get(con.name)
                    if isinstance(reason, str) and reason.strip():
                        used.add(con.name)
                        exempted.append({**finding, "reason": reason.strip()})
                    else:
                        findings.append(finding)
        for name, reason in exceptions.items():
            base = {"module": path.name, "where": "P5_READ1_JUNCTION_EXCEPTIONS",
                    "construct": name, "strand": "-", "segments": "-"}
            if not (isinstance(reason, str) and reason.strip()):
                findings.append({**base, "problem": "exception declared without a reason"})
            elif name not in used:
                findings.append({**base, "problem": "exception matches no junction finding"})
    return {"directory": directory.name, "findings": findings, "skipped": skipped,
            "exempted": exempted}


def protocol_directories() -> list[Path]:
    return sorted(p.parent.parent for p in ROOT.glob("*/tools/*.py")
                  if p.parent.parent.name not in {"tools", "lib", "catalogue"}
                  and p.stem not in SKIP_STEMS
                  and p.parent.parent.parent == ROOT)


def _run(directory: Path, timeout: int) -> dict:
    try:
        done = subprocess.run([sys.executable, __file__, "--worker", str(directory)],
                              capture_output=True, text=True, timeout=timeout, cwd=ROOT)
    except subprocess.TimeoutExpired:
        return {"directory": directory.name, "findings": [], "skipped": ["timed out"]}
    for line in reversed(done.stdout.splitlines()):
        if line.startswith("{"):
            return json.loads(line)
    return {"directory": directory.name, "findings": [],
            "skipped": [f"worker failed: {done.stderr.strip()[-200:]}"]}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("dirs", nargs="*", help="protocol directories (default: all)")
    ap.add_argument("--worker", help=argparse.SUPPRESS)
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--timeout", type=int, default=300, help="seconds per directory")
    ap.add_argument("--show-skipped", action="store_true")
    args = ap.parse_args()
    if args.worker:
        print(json.dumps(scan_directory(Path(args.worker))))
        return 0
    dirs = ([Path(d).resolve() for d in args.dirs] if args.dirs
            else sorted(set(protocol_directories())))
    with cf.ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        results = list(pool.map(lambda d: _run(d, args.timeout), dirs))
    found = [(r["directory"], f) for r in results for f in r["findings"]]
    skipped = [(r["directory"], s) for r in results for s in r["skipped"]]
    exempted = [(r["directory"], f) for r in results for f in r.get("exempted", ())]
    for directory, f in found:
        print(f"{directory}/tools/{f['module']}: {f['construct']!r} ({f['where']}), "
              f"{f['strand']} strand: {f['problem']} at {f['segments']}")
    for directory, f in exempted:
        print(f"exempted {directory}/tools/{f['module']}: {f['construct']!r} ({f['where']}), "
              f"{f['strand']} strand: {f['problem']} at {f['segments']} -- {f['reason']}")
    if args.show_skipped:
        for directory, s in skipped:
            print(f"SKIP {directory}: {s}")
    affected = len({d for d, _ in found})
    print(f"{len(found)} P5/Read 1 junction finding(s) in {affected} director(ies); "
          f"{len(exempted)} exempted; {len(dirs)} scanned, {len(skipped)} module(s) skipped")
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main())
