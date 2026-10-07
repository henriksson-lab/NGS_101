#!/usr/bin/env python3
"""Is the GC effect stronger on genomic DNA than on plasmid? A paired, per-guide test.

Compares two samples of the same library guide-by-guide, so everything guide-intrinsic
cancels -- including the GC-dependent guide-read error rate, which is a property of the
sequence and is therefore identical in both samples.

    delta_g = log2(reads/UMI)_sample2 - log2(reads/UMI)_sample1

regressed on GC. A negative slope means sample 2 loses more reads per molecule at high GC
than sample 1 does.

THREE THINGS THIS HAS TO CONTROL FOR, in order of how badly they bite:

  1. DEPTH. reads-per-UMI grows with sequencing depth, and the growth is not GC-neutral
     because occupancy of the UMI space is abundance-dependent. Comparing an unmatched pair
     nearly doubled the apparent effect here. Match depth with extract.py --max-reads; this
     script warns if the totals differ by more than 5%.
  2. SELECTION. Guide abundance differs between plasmid and cells because cells have been
     under selection. The statistic is reads *per molecule*, so a uniform abundance change
     cancels -- but the estimators are abundance-dependent, so it is controlled explicitly
     here by adding the per-guide abundance shift as a covariate and by stratifying on it.
  3. RUN. The two samples may come from different flowcells. Check with error_control.py
     that the error-vs-GC slope matches before trusting a between-sample difference.

Usage
    python3 paired.py $CHEM_DATA/schmierer/input.matched.guide.tsv \
                      $CHEM_DATA/schmierer/d4.pub.guide.tsv
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
REF = 1


def ols_multi(X: list[list[float]], y: list[float]) -> tuple[list[float], list[float]]:
    """Least squares with an intercept, plus standard errors. Stdlib only."""
    n, k = len(y), len(X[0])
    Xa = [[1.0] + list(x) for x in X]
    p = k + 1
    XtX = [[sum(Xa[i][a] * Xa[i][b] for i in range(n)) for b in range(p)] for a in range(p)]
    Xty = [sum(Xa[i][a] * y[i] for i in range(n)) for a in range(p)]
    aug = [XtX[i][:] + [1.0 if i == j else 0.0 for j in range(p)] + [Xty[i]]
           for i in range(p)]
    for c in range(p):
        piv = max(range(c, p), key=lambda r: abs(aug[r][c]))
        aug[c], aug[piv] = aug[piv], aug[c]
        d = aug[c][c]
        for cc in range(2 * p + 1):
            aug[c][cc] /= d
        for r in range(p):
            if r == c:
                continue
            f = aug[r][c]
            for cc in range(2 * p + 1):
                aug[r][cc] -= f * aug[c][cc]
    beta = [aug[i][2 * p] for i in range(p)]
    resid = [y[i] - sum(beta[a] * Xa[i][a] for a in range(p)) for i in range(n)]
    s2 = sum(e * e for e in resid) / (n - p)
    se = [math.sqrt(s2 * aug[i][p + i]) for i in range(p)]
    return beta, se


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("sample1", type=Path, help="baseline, e.g. the plasmid input")
    ap.add_argument("sample2", type=Path, help="comparison, e.g. the day-4 genomic sample")
    ap.add_argument("--min-distinct", type=int, default=20)
    a = ap.parse_args()

    A = {r["guide"]: r for r in read_guides(a.sample1)}
    B = {r["guide"]: r for r in read_guides(a.sample2)}
    ra, rb = sum(r["reads"] for r in A.values()), sum(r["reads"] for r in B.values())
    print(f"sample1 {a.sample1.name}: {len(A):,} guides, {ra:,} reads")
    print(f"sample2 {a.sample2.name}: {len(B):,} guides, {rb:,} reads")
    skew = abs(ra - rb) / max(ra, rb)
    if skew > 0.05:
        print(f"\n  ** WARNING: depth differs by {100 * skew:.1f}%. reads-per-UMI is "
              f"depth-dependent,\n     so this comparison confounds depth with the "
              f"effect. Re-extract with\n     --max-reads to match before trusting the "
              f"magnitude. **")

    rows = []
    for g, x in A.items():
        y = B.get(g)
        if y is None or x["distinct"] < a.min_distinct or y["distinct"] < a.min_distinct:
            continue
        rows.append(dict(
            gc=x["gc"],
            delta=math.log2(y["reads"] / y["distinct"])
            - math.log2(x["reads"] / x["distinct"]),
            dmol=math.log2(y["distinct"] / x["distinct"])))
    med = st.median(r["delta"] for r in rows)
    for r in rows:
        r["delta"] -= med
    print(f"\npaired on {len(rows):,} guides")

    gc = [r["gc"] for r in rows]
    dy = [r["delta"] for r in rows]
    dm = [r["dmol"] for r in rows]
    b1, s1 = ols_multi([[g] for g in gc], dy)
    b2, s2 = ols_multi([[g, m] for g, m in zip(gc, dm)], dy)
    print(f"\n  delta ~ GC                     GC {b1[1]:+.4f} +/- {s1[1]:.4f}  "
          f"(spearman {spearman(gc, dy):+.4f})")
    print(f"  delta ~ GC + abundance shift   GC {b2[1]:+.4f} +/- {s2[1]:.4f}   "
          f"shift {b2[2]:+.4f} +/- {s2[2]:.4f}")
    print("  (the abundance shift is the selection signal; if the GC term survives it,"
          "\n   the difference is not simply guides dropping out of the cell population)")

    print(f"\npaired difference by GC band, relative to the "
          f"{BANDS[REF][0]:.2f}-{BANDS[REF][1]:.2f} band:")
    cells = []
    for lo, hi in BANDS:
        sub = [r["delta"] for r in rows if lo <= r["gc"] < hi]
        cells.append(st.median(sub) if len(sub) >= 30 else None)
    ref = cells[REF]
    print(f"  {'GC band':12s} {'n':>7s} {'delta':>10s} {'rel. to mid':>12s}")
    for (lo, hi), v in zip(BANDS, cells):
        n = sum(1 for r in rows if lo <= r["gc"] < hi)
        if v is None or ref is None:
            continue
        print(f"  {lo:.2f}-{hi:.2f}   {n:>7,} {v:>+10.4f} {v - ref:>+12.4f}")

    print("\nstratified by the abundance shift (so selection is held roughly constant):")
    rows.sort(key=lambda r: r["dmol"])
    acc: list[list[float]] = [[] for _ in BANDS]
    print(f"  {'shift stratum':18s} " + " ".join(f"{f'{lo:.2f}-{hi:.2f}':>10s}"
                                                 for lo, hi in BANDS))
    for i in range(5):
        ch = rows[i * len(rows) // 5:(i + 1) * len(rows) // 5]
        cs = []
        for lo, hi in BANDS:
            sub = [r["delta"] for r in ch if lo <= r["gc"] < hi]
            cs.append(st.median(sub) if len(sub) >= 30 else None)
        rf = cs[REF]
        if rf is None:
            continue
        out = []
        for j, c in enumerate(cs):
            if c is None:
                out.append(f"{'-':>10s}")
            else:
                acc[j].append(c - rf)
                out.append(f"{c - rf:>+10.4f}")
        print(f"  {ch[0]['dmol']:+.2f}..{ch[-1]['dmol']:+.2f}      " + " ".join(out))
    print(f"  {'mean':18s} " + " ".join(f"{st.mean(c):>+10.4f}" if c else f"{'-':>10s}"
                                        for c in acc))
    if acc[-1]:
        neg = sum(1 for v in acc[-1] if v < 0)
        print(f"\nhigh-GC band: mean {st.mean(acc[-1]):+.4f}, negative in "
              f"{neg}/{len(acc[-1])} strata of the abundance shift")
    return 0


if __name__ == "__main__":
    sys.exit(main())
