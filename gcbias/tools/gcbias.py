#!/usr/bin/env python3
"""Compare observed read counts to what the lineage-UMI distribution predicts, vs GC%.

Procedure
    1. Fit ONE global (mu, r) across all guides, from the pooled UMI read-count moments.
       Global on purpose: it carries no per-guide information, so the prediction it makes
       for a guide uses only that guide's molecule count.
    2. Per guide, invert the observed distinct-UMI count to a molecule estimate N_hat
       (correcting label collisions and unseen molecules), then predict reads = N_hat * mu.
    3. bias = log2(observed reads / predicted reads). Positive = over-amplified.
    4. Regress bias on GC%. Report the naive log2(reads/distinct) alongside, because the
       difference between the two is the whole argument.

Usage
    python3 gcbias.py $CHEM_DATA/schmierer/input.25000000.guide.tsv --max-occupancy 0.8
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import umimodel as um  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from qualitative import read_guides, spearman  # noqa: E402


def ols(xs: list[float], ys: list[float]) -> tuple[float, float, float, float]:
    """slope, intercept, r, stderr(slope) -- stdlib only."""
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    syy = sum((y - my) ** 2 for y in ys)
    if sxx == 0:
        return float("nan"), my, float("nan"), float("nan")
    b = sxy / sxx
    a = my - b * mx
    r = sxy / math.sqrt(sxx * syy) if syy > 0 else float("nan")
    resid = syy - b * sxy
    se = math.sqrt(resid / (n - 2) / sxx) if n > 2 and resid > 0 else float("nan")
    return b, a, r, se


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("guide_tsv", type=Path)
    ap.add_argument("--umi-len", type=int, default=6)
    ap.add_argument("--max-occupancy", type=float, default=0.90,
                    help="drop guides whose UMI space is fuller than this")
    ap.add_argument("--min-reads", type=int, default=100)
    ap.add_argument("--min-distinct", type=int, default=20)
    ap.add_argument("--guide-umi", type=Path, default=None,
                    help="*.guide_umi.tsv -- enables the non-parametric (Chao1) estimator, "
                         "which is the one to use when the parametric fit is rejected")
    a = ap.parse_args()

    space = 4 ** a.umi_len
    allrows = read_guides(a.guide_tsv)
    rows = [r for r in allrows
            if r["reads"] >= a.min_reads and r["distinct"] >= a.min_distinct
            and r["distinct"] / space <= a.max_occupancy]
    print(f"guides: {len(allrows):,} -> {len(rows):,} after filters "
          f"(reads>={a.min_reads}, distinct>={a.min_distinct}, "
          f"occupancy<={a.max_occupancy})")
    if len(rows) < 20:
        raise SystemExit("too few guides survive the filters to say anything")

    # per-guide singleton/doubleton counts for the non-parametric route
    f12: dict[str, tuple[int, int]] = {}
    if a.guide_umi and a.guide_umi.exists():
        from collections import defaultdict
        c1: dict[str, int] = defaultdict(int)
        c2: dict[str, int] = defaultdict(int)
        with a.guide_umi.open() as fh:
            fh.readline()
            for line in fh:
                g, _, n = line.rstrip("\n").rsplit("\t", 2)
                k = int(n)
                if k == 1:
                    c1[g] += 1
                elif k == 2:
                    c2[g] += 1
        f12 = {g: (c1.get(g, 0), c2.get(g, 0)) for g in {*c1, *c2}}
        print(f"singleton/doubleton counts loaded for {len(f12):,} guides")

    # ---------------------------------------------------------------- 1. global fit
    tot_d = sum(r["distinct"] for r in rows)
    pooled_m1 = sum(r["m1"] * r["distinct"] for r in rows) / tot_d
    pooled_m2 = sum(r["m2"] * r["distinct"] for r in rows) / tot_d
    mean_occ = sum(r["distinct"] for r in rows) / (len(rows) * space)
    fit = um.fit_moments(space, mean_occ * space, pooled_m1, pooled_m2)
    print(f"\nglobal fit on pooled moments (m1={pooled_m1:.3f}, m2={pooled_m2:.1f}, "
          f"mean occupancy={mean_occ:.4f}):")
    if not fit.ok:
        print(f"  FAILED: {fit.reason}")
        print("  Falling back to a Poisson amplification model (r -> inf). The GC")
        print("  comparison still runs, but the dispersion is not identified.")
        mu, r = pooled_m1, float("inf")
    else:
        mu, r = fit.mu, fit.r
        print(f"  mu (reads per molecule) = {mu:.3f}")
        print(f"  r  (amplification shape) = {r:.3f}"
              + ("   [Poisson limit]" if math.isinf(r) else ""))
        print(f"  implied P(molecule unseen) = {um.p0_nb(mu, r):.4f}")
        print(f"  {fit.reason}" if fit.reason else "")

    # ---------------------------------------------------------------- 2/3. per-guide bias
    recs = []
    for g in rows:
        n_hat = um.estimate_molecules(g["distinct"], space, mu, r)
        pred = um.expected_reads(n_hat, mu)
        if not math.isfinite(pred) or pred <= 0:
            continue
        rec = dict(
            guide=g["guide"], gc=g["gc"], reads=g["reads"], distinct=g["distinct"],
            occupancy=g["distinct"] / space, n_hat=n_hat, pred=pred,
            bias=math.log2(g["reads"] / pred),
            naive=math.log2(g["reads"] / g["distinct"]),
            rel_se=um.molecules_rel_se(g["distinct"], space, mu, r),
            np_n_hat=float("nan"), np_bias=float("nan"), coverage=float("nan"))
        if g["guide"] in f12:
            f1, f2 = f12[g["guide"]]
            nn = um.nonparametric_molecules(g["distinct"], f1, f2, space)
            rec["np_n_hat"] = nn
            rec["coverage"] = um.good_turing_coverage(g["reads"], f1)
            if math.isfinite(nn) and nn > 0:
                rec["np_bias"] = math.log2(g["reads"] / nn)   # reads per molecule, log2
        recs.append(rec)
    # Expected reads from the non-parametric molecule estimate, scaled by ONE global
    # reads-per-molecule so the totals match. That makes the identity line the reference:
    # a guide off the line is amplifying differently from the library average, rather than
    # reflecting a global offset in the fitted mu.
    np_ok = [x for x in recs if math.isfinite(x["np_n_hat"]) and x["np_n_hat"] > 0]
    if np_ok:
        mu_bar = sum(x["reads"] for x in np_ok) / sum(x["np_n_hat"] for x in np_ok)
        for x in recs:
            x["np_pred"] = (x["np_n_hat"] * mu_bar
                            if math.isfinite(x["np_n_hat"]) else float("nan"))
        print(f"\nglobal reads per molecule (Chao1 route): mu_bar = {mu_bar:.4f}")
    else:
        for x in recs:
            x["np_pred"] = float("nan")

    print(f"\nmodelled {len(recs):,} guides")
    if f12:
        cov = sorted(x["coverage"] for x in recs if math.isfinite(x["coverage"]))
        nn = [x["np_n_hat"] for x in recs if math.isfinite(x["np_n_hat"])]
        if cov:
            print(f"Good-Turing coverage: median={cov[len(cov) // 2]:.3f} "
                  f"(low = too shallow for any dedup claim)")
        if nn:
            srt = sorted(nn)
            print(f"Chao1+collision molecules/guide: median={srt[len(srt) // 2]:,.0f} "
                  f"(paper: ~3,355 lineages/guide)")

    # centre both statistics so slopes, not offsets, are compared
    for key in ("bias", "naive", "np_bias"):
        vals = sorted(x[key] for x in recs if math.isfinite(x[key]))
        if not vals:
            continue
        med = vals[len(vals) // 2]
        for x in recs:
            x[key] -= med

    gc = [x["gc"] for x in recs]
    print("\n" + "=" * 68)
    print("GC% REGRESSION   (slope is in log2 units per unit GC fraction)")
    print("=" * 68)
    print(f"{'statistic':26s} {'slope':>9s} {'per 10%GC':>10s} {'r':>7s} {'se':>8s} "
          f"{'spearman':>9s}")
    for key, name in (("naive", "naive log2(reads/distinct)"),
                      ("bias", "NB model log2(obs/pred)"),
                      ("np_bias", "Chao1 log2(reads/molecule)")):
        pairs = [(x["gc"], x[key]) for x in recs if math.isfinite(x[key])]
        if len(pairs) < 20:
            continue
        xs = [p[0] for p in pairs]
        ys = [p[1] for p in pairs]
        b, _, rr, se = ols(xs, ys)
        print(f"{name:26s} {b:>9.4f} {b / 10:>10.4f} {rr:>7.3f} {se:>8.4f} "
              f"{spearman(xs, ys):>9.4f}")
    print("\nIf these two disagree, the naive one is wrong: it reads the UMI-space")
    print("saturation curve as a GC effect. Check occupancy spread below.")

    occs = sorted(x["occupancy"] for x in recs)
    print(f"\noccupancy of modelled guides: min={occs[0]:.3f} "
          f"median={occs[len(occs) // 2]:.3f} max={occs[-1]:.3f}")
    ses = sorted(x["rel_se"] for x in recs if math.isfinite(x["rel_se"]))
    if ses:
        print(f"relative SE of N_hat:         median={ses[len(ses) // 2]:.3f} "
              f"p95={ses[int(0.95 * len(ses))]:.3f}")

    stat = "np_bias" if f12 else "bias"
    label = "Chao1" if f12 else "NB model"
    print(f"\nbias by GC quintile ({label} statistic):")
    srt = sorted((x for x in recs if math.isfinite(x[stat])), key=lambda x: x["gc"])
    for i in range(5):
        ch = srt[i * len(srt) // 5:(i + 1) * len(srt) // 5]
        if not ch:
            continue
        bs = sorted(x[stat] for x in ch)
        print(f"  Q{i + 1} GC {ch[0]['gc']:.2f}-{ch[-1]['gc']:.2f}  n={len(ch):>6,}  "
              f"median bias={bs[len(bs) // 2]:+.4f}  "
              f"IQR=[{bs[len(bs) // 4]:+.3f},{bs[3 * len(bs) // 4]:+.3f}]")

    out = a.guide_tsv.with_suffix("").with_suffix("")
    p = Path(f"{out}.gcbias.tsv")
    with p.open("w") as fh:
        fh.write("guide\tgc\treads\tdistinct_umi\toccupancy\tn_hat\tpred_reads\t"
                 "bias_log2\tnaive_log2\trel_se\tnp_n_hat\tnp_pred_reads\t"
                 "np_bias_log2\tcoverage\n")
        for x in sorted(recs, key=lambda y: y["gc"]):
            fh.write(f"{x['guide']}\t{x['gc']:.4f}\t{x['reads']}\t{x['distinct']}\t"
                     f"{x['occupancy']:.6f}\t{x['n_hat']:.2f}\t{x['pred']:.2f}\t"
                     f"{x['bias']:.6f}\t{x['naive']:.6f}\t{x['rel_se']:.6f}\t"
                     f"{x['np_n_hat']:.2f}\t{x['np_pred']:.2f}\t"
                     f"{x['np_bias']:.6f}\t{x['coverage']:.6f}\n")
    print(f"\nwrote {p}")
    with open(f"{out}.globalfit.tsv", "w") as fh:
        fh.write("param\tvalue\n")
        fh.write(f"umi_len\t{a.umi_len}\nspace\t{space}\nmu\t{mu}\nr\t{r}\n"
                 f"pooled_m1\t{pooled_m1}\npooled_m2\t{pooled_m2}\n"
                 f"n_guides\t{len(recs)}\n")
    print(f"wrote {out}.globalfit.tsv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
