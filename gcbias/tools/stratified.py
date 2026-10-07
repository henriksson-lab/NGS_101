#!/usr/bin/env python3
"""Model-free GC test: hold the molecule count fixed, ask whether GC changes the reads.

Why this exists. Guide abundance is non-uniform in the plasmid library, and every estimator
in gcbias.py (Chao1, the NB inversion, the collision correction) is a function of abundance.
So a GC effect could in principle be an abundance effect laundered through the estimator.

This test uses no model at all. Guides are stratified into narrow bins of *observed distinct
UMI count* -- a direct molecule proxy -- and within each stratum log2(reads) is regressed on
GC. If GC-rich guides with the same number of distinct UMIs get fewer reads, that is
amplification, not an artefact of molecule estimation.

It is also CONSERVATIVE. If GC-rich guides really do amplify less, they lose reads *and*
lose detected UMIs, so at fixed distinct-UMI count they carry slightly more true molecules
and should show more reads, not fewer. Any deficit measured here is therefore a lower bound
on the real effect.

What it does NOT separate: the UMI labels a cloning event, so "reads per molecule" spans
bacterial propagation, plasmid copy number, miniprep and PCR. See the README.

Usage
    python3 stratified.py $CHEM_DATA/schmierer/input.collapsed.guide.tsv
"""
from __future__ import annotations

import argparse
import math
import statistics as st
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "lib"))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

from gcbias import ols  # noqa: E402
from qualitative import read_guides, spearman  # noqa: E402

BANDS = [(0.0, 0.45), (0.45, 0.55), (0.55, 0.60), (0.60, 0.70), (0.70, 1.01)]
REF = 1  # index of the reference band, 0.45-0.55


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("guide_tsv", type=Path)
    ap.add_argument("--strata", type=int, default=10)
    ap.add_argument("--min-per-cell", type=int, default=30)
    a = ap.parse_args()

    rows = read_guides(a.guide_tsv)
    rows.sort(key=lambda r: r["distinct"])
    n = len(rows)
    print(f"{n:,} guides, {a.strata} strata of observed distinct-UMI count\n")

    print("A. slope of log2(reads) on GC, within each stratum")
    print(f"   {'stratum (distinct UMI)':24s} {'n':>6s} {'med reads':>10s} "
          f"{'slope':>18s} {'spearman':>9s}")
    slopes = []
    for i in range(a.strata):
        ch = rows[i * n // a.strata:(i + 1) * n // a.strata]
        if len(ch) < 200:
            continue
        gc = [x["gc"] for x in ch]
        y = [math.log2(x["reads"]) for x in ch]
        b, _, _, se = ols(gc, y)
        slopes.append(b)
        lab = f"{ch[0]['distinct']}-{ch[-1]['distinct']}"
        print(f"   {lab:24s} {len(ch):>6,} "
              f"{st.median([x['reads'] for x in ch]):>10.0f} "
              f"{b:>+12.4f}±{se:.4f} {spearman(gc, y):>+9.4f}")
    if slopes:
        neg = sum(1 for s in slopes if s < 0)
        print(f"\n   mean slope {st.mean(slopes):+.4f} log2 per unit GC; "
              f"{neg}/{len(slopes)} strata negative")

    print("\nB. by GC band within strata, relative to the "
          f"{BANDS[REF][0]:.2f}-{BANDS[REF][1]:.2f} band")
    print(f"   {'stratum':24s} " + " ".join(f"{f'{lo:.2f}-{hi:.2f}':>10s}"
                                            for lo, hi in BANDS))
    col_acc: list[list[float]] = [[] for _ in BANDS]
    step = max(1, a.strata // 5)
    for i in range(0, a.strata, step):
        ch = rows[i * n // a.strata:min(a.strata, i + step) * n // a.strata]
        cells = []
        for lo, hi in BANDS:
            sub = [x for x in ch if lo <= x["gc"] < hi]
            cells.append(math.log2(st.median([x["reads"] for x in sub]))
                         if len(sub) >= a.min_per_cell else None)
        ref = cells[REF]
        if ref is None:
            continue
        out = []
        for j, c in enumerate(cells):
            if c is None:
                out.append(f"{'-':>10s}")
            else:
                col_acc[j].append(c - ref)
                out.append(f"{c - ref:>+10.4f}")
        lab = f"{ch[0]['distinct']}-{ch[-1]['distinct']}"
        print(f"   {lab:24s} " + " ".join(out))

    print(f"\n   {'mean over strata':24s} " +
          " ".join(f"{st.mean(c):>+10.4f}" if c else f"{'-':>10s}" for c in col_acc))
    print("\n   (log2 reads relative to the mid-GC band at the same molecule count;"
          "\n    negative = fewer reads than mid-GC guides with the same molecule count)")

    hi_band = col_acc[-1]
    lo_band = col_acc[0]
    print("\nC. verdict")
    if hi_band:
        consistent = sum(1 for v in hi_band if v < 0)
        print(f"   high GC (>{BANDS[-1][0]:.2f}): mean {st.mean(hi_band):+.4f}, "
              f"negative in {consistent}/{len(hi_band)} strata")
    if lo_band:
        consistent = sum(1 for v in lo_band if v < 0)
        print(f"   low  GC (<{BANDS[0][1]:.2f}): mean {st.mean(lo_band):+.4f}, "
              f"negative in {consistent}/{len(lo_band)} strata")
    print("   A deficit that holds across strata is not an abundance artefact.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
