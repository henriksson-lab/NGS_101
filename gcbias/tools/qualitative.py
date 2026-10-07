#!/usr/bin/env python3
"""Qualitative diagnostics -- run this BEFORE fitting anything.

The point of looking first is that every quantitative statement below depends on
assumptions that are visible by eye and invisible in a p-value: whether the UMI space is
saturated, whether error reads have filled it, whether GC correlates with abundance (which
would make a GC "bias" out of pure biology), and whether any effect lives in the extremes
rather than the bulk.

Usage
    python3 qualitative.py $CHEM_DATA/schmierer/input.25000000.guide.tsv \
        --guide-umi $CHEM_DATA/schmierer/input.25000000.guide_umi.tsv

Prints ASCII, and writes TSVs for plotting in R (see ../R/).
"""
from __future__ import annotations

import argparse
import math
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import umimodel as um  # noqa: E402


def read_guides(p: Path) -> list[dict]:
    rows = []
    with p.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            v = line.rstrip("\n").split("\t")
            d = dict(zip(hdr, v))
            rows.append(dict(guide=d["guide"], gc=float(d["gc"]), reads=int(d["reads"]),
                             distinct=int(d["distinct_umi"]), m1=float(d["m1"]),
                             m2=float(d["m2"])))
    return rows


def bar(frac: float, width: int = 40) -> str:
    n = max(0, min(width, int(round(frac * width))))
    return "#" * n + "." * (width - n)


def hist(values: list[float], bins: int, lo: float | None = None,
         hi: float | None = None, label: str = "", fmt: str = "{:.2f}") -> None:
    if not values:
        print("  (no data)")
        return
    lo = min(values) if lo is None else lo
    hi = max(values) if hi is None else hi
    if hi <= lo:
        hi = lo + 1e-9
    counts = [0] * bins
    for v in values:
        k = min(bins - 1, max(0, int((v - lo) / (hi - lo) * bins)))
        counts[k] += 1
    top = max(counts) or 1
    for i, c in enumerate(counts):
        a = lo + (hi - lo) * i / bins
        b = lo + (hi - lo) * (i + 1) / bins
        print(f"  {fmt.format(a)}..{fmt.format(b)} {bar(c / top, 34)} {c:>7,}")


def quantiles(xs: list[float], qs=(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99)) -> dict:
    s = sorted(xs)
    return {q: s[min(len(s) - 1, int(q * len(s)))] for q in qs}


def spearman(xs: list[float], ys: list[float]) -> float:
    """Rank correlation, no scipy needed."""
    def ranks(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = ranks(xs), ranks(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return num / den if den else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("guide_tsv", type=Path)
    ap.add_argument("--guide-umi", type=Path, default=None)
    ap.add_argument("--umi-len", type=int, default=6)
    ap.add_argument("--per-base-error", type=float, default=0.005)
    ap.add_argument("--min-reads", type=int, default=0)
    a = ap.parse_args()

    space = 4 ** a.umi_len
    rows = [r for r in read_guides(a.guide_tsv) if r["reads"] >= a.min_reads]
    if not rows:
        raise SystemExit("no guides pass --min-reads")
    out = a.guide_tsv.with_suffix("")

    print("=" * 72)
    print(f"QUALITATIVE DIAGNOSTICS   {a.guide_tsv.name}")
    print(f"{len(rows):,} guides, UMI length {a.umi_len} -> label space {space:,}")
    print("=" * 72)

    # ---------------------------------------------------------------- 1. is it saturated?
    print("\n[1] UMI SPACE OCCUPANCY per guide   (distinct UMIs / label space)")
    print("    The model's inversion needs headroom. Near 1.0, dedup is meaningless.")
    occ = [r["distinct"] / space for r in rows]
    hist(occ, 12, 0.0, 1.0)
    q = quantiles(occ)
    print("    quantiles: " + "  ".join(f"p{int(k * 100)}={v:.3f}" for k, v in q.items()))
    n_sat = sum(1 for o in occ if o > 0.9)
    print(f"    guides above 0.90 occupancy: {n_sat:,} ({100 * n_sat / len(occ):.1f}%)"
          + ("   <-- these cannot be deduplicated" if n_sat else ""))

    # ------------------------------------------------- 2. has error filled the UMI space?
    print(f"\n[2] IS OCCUPANCY REAL, OR SEQUENCING ERROR?  "
          f"(assuming {a.per_base_error:.1%} per-base)")
    thresh = um.reads_before_error_fills_space(space, a.umi_len, a.per_base_error)
    print(f"    P(UMI carries an error)          {um.umi_error_rate(a.umi_len, a.per_base_error):.4f}")
    print(f"    reads/guide at which error alone fills the space: {thresh:,.0f}")
    over = sum(1 for r in rows if r["reads"] > thresh)
    print(f"    guides above that depth: {over:,} ({100 * over / len(rows):.1f}%)")
    print("    -> above it, occupancy measures depth, not molecules. Collapse UMIs"
          "\n       (extract.py --collapse 1) and/or subsample, then re-check.")

    # ------------------------------------------- 3. the confound: does GC track abundance?
    print("\n[3] THE CONFOUND: does GC% correlate with raw abundance?")
    print("    If it does, any reads-vs-UMI trend may be biology, not amplification.")
    gc = [r["gc"] for r in rows]
    reads = [math.log2(r["reads"]) for r in rows]
    dis = [r["distinct"] for r in rows]
    print(f"    Spearman(GC, log2 reads)        {spearman(gc, reads):+.4f}")
    print(f"    Spearman(GC, distinct UMIs)     {spearman(gc, [float(x) for x in dis]):+.4f}")
    print(f"    Spearman(GC, reads per UMI)     "
          f"{spearman(gc, [r['reads'] / r['distinct'] for r in rows]):+.4f}")

    # ------------------------------------------------ 4. the duplication spectrum by GC
    print("\n[4] READS PER UMI, BY GC QUINTILE   (the raw, uncorrected ratio)")
    print("    A clean amplification effect shows as a monotone shift across quintiles.")
    srt = sorted(rows, key=lambda r: r["gc"])
    nq = 5
    print(f"    {'quintile':10s} {'GC range':14s} {'guides':>7s} {'reads/UMI':>10s} "
          f"{'occupancy':>10s} {'median reads':>13s}")
    qstats = []
    for i in range(nq):
        chunk = srt[i * len(srt) // nq:(i + 1) * len(srt) // nq]
        if not chunk:
            continue
        rpu = sorted(c["reads"] / c["distinct"] for c in chunk)
        med_rpu = rpu[len(rpu) // 2]
        oc = sorted(c["distinct"] / space for c in chunk)
        rd = sorted(c["reads"] for c in chunk)
        print(f"    Q{i + 1:<9d} {chunk[0]['gc']:.2f}-{chunk[-1]['gc']:.2f}     "
              f"{len(chunk):>7,} {med_rpu:>10.3f} {oc[len(oc) // 2]:>10.3f} "
              f"{rd[len(rd) // 2]:>13,}")
        qstats.append(dict(q=i + 1, gc_lo=chunk[0]["gc"], gc_hi=chunk[-1]["gc"],
                           n=len(chunk), reads_per_umi=med_rpu,
                           occupancy=oc[len(oc) // 2], median_reads=rd[len(rd) // 2]))

    # --------------------------------------------------------------- 5. the extreme tails
    print("\n[5] THE EXTREMES   (where the effect was reported to show up)")
    print("    Bulk GC may look flat while the tails move. Compare ends, not slope.")
    for lo, hi, name in ((0.0, 0.30, "GC < 30%"), (0.30, 0.70, "GC 30-70% (bulk)"),
                         (0.70, 1.01, "GC > 70%")):
        sub = [r for r in rows if lo <= r["gc"] < hi]
        if not sub:
            print(f"    {name:18s} (none)")
            continue
        rpu = sorted(r["reads"] / r["distinct"] for r in sub)
        print(f"    {name:18s} n={len(sub):>6,}  median reads/UMI={rpu[len(rpu) // 2]:7.3f}"
              f"  median reads={sorted(r['reads'] for r in sub)[len(sub) // 2]:>8,}")

    # ------------------------------------------------------- 6. per-UMI read distribution
    if a.guide_umi and a.guide_umi.exists():
        print("\n[6] PER-UMI READ-COUNT SPECTRUM   (pooled; the shape the model fits)")
        spec = Counter()
        with a.guide_umi.open() as fh:
            fh.readline()
            for line in fh:
                spec[int(line.rsplit("\t", 1)[1])] += 1
        tot = sum(spec.values())
        print(f"    {tot:,} (guide, UMI) pairs")
        top = max(spec.values())
        for k in range(1, min(16, max(spec) + 1)):
            c = spec.get(k, 0)
            print(f"      {k:>3d} read(s) {bar(c / top, 34)} {c:>9,} "
                  f"({100 * c / tot:5.2f}%)")
        tail = sum(v for k, v in spec.items() if k >= 16)
        print(f"      >=16 reads {bar(tail / top, 34)} {tail:>9,} "
              f"({100 * tail / tot:5.2f}%)")
        with open(f"{out}.umi_spectrum.tsv", "w") as fh:
            fh.write("reads_per_umi\tn_umis\n")
            for k in sorted(spec):
                fh.write(f"{k}\t{spec[k]}\n")
        print(f"    wrote {out}.umi_spectrum.tsv")

    with open(f"{out}.gc_quintiles.tsv", "w") as fh:
        fh.write("quintile\tgc_lo\tgc_hi\tn\treads_per_umi\toccupancy\tmedian_reads\n")
        for s in qstats:
            fh.write(f"{s['q']}\t{s['gc_lo']:.4f}\t{s['gc_hi']:.4f}\t{s['n']}\t"
                     f"{s['reads_per_umi']:.6f}\t{s['occupancy']:.6f}\t{s['median_reads']}\n")
    print(f"\nwrote {out}.gc_quintiles.tsv")
    print("\nRead [1] and [2] before trusting anything in [3]-[5].")
    return 0


if __name__ == "__main__":
    sys.exit(main())
