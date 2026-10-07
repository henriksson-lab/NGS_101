#!/usr/bin/env python3
"""Show what the naive reads/UMI ratio does when there is NO GC effect at all.

This is the calibration run that justifies the whole model. We simulate a library with:

  * no GC-dependent amplification whatsoever -- every guide shares one (mu, r);
  * abundance that correlates with GC, which is true in real libraries (GC-rich guides
    clone and grow differently, and GC correlates with gene class).

Under those conditions the *truth* is a zero GC slope. The naive statistic nevertheless
reports a large one, purely because distinct-UMI counts saturate at the label space and
abundance is confounded with GC. The model statistic should recover ~0.

Then the reverse: inject a real GC-dependent amplification effect and check the model
recovers its size while the naive statistic mis-states it.

Usage
    python3 simulate.py                      # 6-nt RSL, Schmierer-like
    python3 simulate.py --umi-len 10         # Michlits-like
"""
from __future__ import annotations

import argparse
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "lib"))
import umimodel as um  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gcbias import ols  # noqa: E402


def run(n_guides: int, space: int, mu: float, r: float, gc_amp_log2: float,
        abundance_gc_coupling: float, n_mid: int, rng: random.Random) -> dict:
    """One synthetic library. gc_amp_log2 is the TRUE log2 amplification change per unit GC."""
    gcs, naive, model, occs = [], [], [], []
    for _ in range(n_guides):
        gc = min(0.95, max(0.05, rng.gauss(0.5, 0.12)))
        # abundance deliberately confounded with GC, but amplification maybe not
        n_mol = max(20, int(n_mid * math.exp(abundance_gc_coupling * (gc - 0.5))))
        mu_g = mu * (2.0 ** (gc_amp_log2 * (gc - 0.5)))
        d = um.simulate_guide(n_mol, space, mu_g, r, rng)
        if len(d) < 10:
            continue
        D, R = len(d), sum(d.values())
        if D >= space:
            continue
        gcs.append(gc)
        naive.append(math.log2(R / D))
        model.append(math.log2(R / um.predict_reads(D, space, mu, r)))
        occs.append(D / space)
    nb, _, _, nse = ols(gcs, naive)
    mb, _, _, mse = ols(gcs, model)
    return dict(n=len(gcs), naive_slope=nb, naive_se=nse, model_slope=mb, model_se=mse,
                occ_med=sorted(occs)[len(occs) // 2], occ_max=max(occs),
                truth=gc_amp_log2)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--umi-len", type=int, default=6)
    ap.add_argument("--guides", type=int, default=600)
    ap.add_argument("--mu", type=float, default=20.0)
    ap.add_argument("--r", type=float, default=5.0)
    ap.add_argument("--seed", type=int, default=11)
    a = ap.parse_args()
    space = 4 ** a.umi_len
    rng = random.Random(a.seed)

    print("=" * 78)
    print(f"CALIBRATION   UMI length {a.umi_len} -> label space {space:,};  "
          f"mu={a.mu}, r={a.r}")
    print("=" * 78)
    print("\nA. NO true GC effect on amplification, but abundance correlates with GC.")
    print("   Truth = 0.0. Any nonzero slope is an artefact of the statistic.\n")
    print(f"   {'lineages/guide':>15s} {'med occ':>8s} {'naive slope':>13s} "
          f"{'model slope':>13s}")
    for n_mid in (200, 800, 2000, 4000, 12000):
        res = run(a.guides, space, a.mu, a.r, 0.0, 1.2, n_mid, rng)
        if res["n"] < 50:
            print(f"   {n_mid:>15,} -- too many guides saturate the space entirely")
            continue
        print(f"   {n_mid:>15,} {res['occ_med']:>8.3f} "
              f"{res['naive_slope']:>+9.3f}±{res['naive_se']:.3f} "
              f"{res['model_slope']:>+9.3f}±{res['model_se']:.3f}")
    print("\n   The naive slope grows with occupancy even though the truth stays 0.")
    print("   That is the false GC signal saturation manufactures.")

    print("\nB. A REAL GC effect injected (true slope -1.0 log2 per unit GC,")
    print("   i.e. GC-rich guides amplify less). Can each statistic recover it?\n")
    print(f"   {'lineages/guide':>15s} {'med occ':>8s} {'naive slope':>13s} "
          f"{'model slope':>13s} {'truth':>7s}")
    for n_mid in (200, 800, 2000, 4000):
        res = run(a.guides, space, a.mu, a.r, -1.0, 1.2, n_mid, rng)
        if res["n"] < 50:
            continue
        print(f"   {n_mid:>15,} {res['occ_med']:>8.3f} "
              f"{res['naive_slope']:>+9.3f}±{res['naive_se']:.3f} "
              f"{res['model_slope']:>+9.3f}±{res['model_se']:.3f} {-1.0:>7.1f}")

    print("\nC. Estimator recovery of molecule counts, same model.\n")
    print(f"   {'N true':>9s} {'occupancy':>10s} {'N_hat':>10s} {'error':>8s}")
    for n_mol in (100, 500, 1500, 3000, 6000, 12000):
        d = um.simulate_guide(n_mol, space, a.mu, a.r, rng)
        D = len(d)
        if D >= space:
            print(f"   {n_mol:>9,} {D / space:>10.3f}    space fully occupied -- "
                  f"no estimate possible")
            continue
        nh = um.estimate_molecules(D, space, a.mu, a.r)
        print(f"   {n_mol:>9,} {D / space:>10.3f} {nh:>10,.0f} "
              f"{100 * (nh - n_mol) / n_mol:>+7.1f}%")

    print("\nTakeaway: report the model statistic, and report the occupancy range it was")
    print("computed over. Quote the naive one only to show the difference.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
