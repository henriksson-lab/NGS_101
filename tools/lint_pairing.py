#!/usr/bin/env python3
"""
Lint every generated page for duplex drawings whose columns do not base-pair.

Two consecutive lines inside a <pre> block, one drawn 5'->3' and the next 3'->5' (read off
their 5'-/-3' and 3'-/-5' end labels), are taken to be a duplex: every column where both
carry a base must satisfy chemdraw.bases_pair. That catches the hand-indent class of error
(an oligo-dT over the cDNA, a GGG nowhere near its CCC) whichever chemistry it is in.

New drawings should be built with chemdraw.Scene, which cannot produce these at all.

Run:  python3 tools/lint_pairing.py [files...]     (default: every *.html under the repo)
Exit status 1 if anything is found.
"""
from __future__ import annotations

import html
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))
from chemdraw import bases_pair  # noqa: E402

SKIP = {"ref", "to_debug", "_data", "gcbias"}          # archived exemplars, not ours
SEQ_TOKEN = re.compile(r"[ACGTUNWSRYKMBDHVXIacgtunwsrykmbdhvxi.]{3,}")
WORD = re.compile(r"[A-Za-z.]+")


def orientation(line: str) -> str | None:
    fwd = "5'-" in line or "-3'" in line
    rev = "3'-" in line or "-5'" in line
    if fwd == rev:
        return None
    return "53" if fwd else "35"


NOT_BASES = re.compile(r"/[^/\s]+/|\[[^\]]*\]|\([^)]*\)")   # /5Phos/, [8bp UMI], (A)n


def cells(line: str) -> dict[int, str]:
    out = {}
    line = NOT_BASES.sub(lambda m: " " * len(m.group()), line)
    # only what lies between the strand's own end labels is sequence: a gutter label
    # ("cDNA") before 5'-/3'- or a note after -3'/-5' is not
    heads = [i for i in (line.find("5'-"), line.find("3'-")) if i >= 0]
    if heads:
        h = min(heads) + 3
        line = " " * h + line[h:]
    tails = [i for i in (line.rfind("-3'"), line.rfind("-5'")) if i >= 0]
    if tails:
        line = line[:max(tails)]
    for m in WORD.finditer(line):
        w = m.group()
        if SEQ_TOKEN.fullmatch(w) and not set(w) <= {"."}:
            for i, ch in enumerate(m.group()):
                out[m.start() + i] = ch
    return out


def lint_text(text: str) -> list[tuple[str, str, str, int]]:
    """-> [(caption, upper line, lower line, n bad columns)]"""
    found = []
    for block in re.findall(r"<pre>(.*?)</pre>", text, flags=re.S):
        cap = re.search(r"<i>(.*?)</i>", block, flags=re.S)
        caption = html.unescape(re.sub(r"<[^>]+>", "", cap.group(1)))[:90] if cap else ""
        # <unp> marks segments a Scene drew deliberately unpaired (a 5' flap): blank them
        block = re.sub(r"<unp>(.*?)</unp>", lambda m: " " * len(
            html.unescape(re.sub(r"<[^>]+>", "", m.group(1)))), block, flags=re.S)
        plain = html.unescape(re.sub(r"<[^>]+>", "", block))
        lines = plain.split("\n")
        for a, b in zip(lines, lines[1:]):
            oa, ob = orientation(a), orientation(b)
            if not oa or not ob or oa == ob:
                continue
            ca, cb = cells(a), cells(b)
            both = set(ca) & set(cb)
            bad = [c for c in both if not bases_pair(ca[c], cb[c])]
            if bad:
                marks = "".join("^" if i in bad else " " for i in range(max(bad) + 1))
                found.append((caption, a.rstrip(), b.rstrip() + "\n    " + marks, len(bad)))
    return found


def main(argv: list[str]) -> int:
    files = [Path(a).resolve() for a in argv] or sorted(
        p for p in ROOT.rglob("*.html") if not (set(p.relative_to(ROOT).parts) & SKIP))
    total = 0
    for f in files:
        for caption, a, b, n in lint_text(f.read_text()):
            total += 1
            print(f"{f.relative_to(ROOT)}: {n} unpaired column(s)  [{caption}]\n    {a}\n    {b}\n")
    print(f"{total} problem duplex(es)")
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
