#!/usr/bin/env python3
"""Thin a guide-UMI table to a target reads-per-UMI, so samples can be compared fairly.

WHY THIS EXISTS -- the single most important methodological point in this directory.

The working statistic is log2(reads) at fixed observed distinct UMIs, i.e. essentially
log(R/D). But R/D is not the per-molecule amplification rate mu. Molecules that received
zero reads are invisible, so what we observe is the ZERO-TRUNCATED mean

    T(mu) = mu / (1 - p0(mu))

which is monotone in mu but compressed, and the compression depends on where the sample
sits:

    dlogT/dlogmu  ->  0   as mu -> 0   (shallow: a real effect is squashed to nothing)
    dlogT/dlogmu  ->  1   as mu -> inf (deep: the effect shows in full)

So two samples with the SAME underlying GC effect report different effect sizes purely
because they sit at different reads-per-molecule. Matching total read depth does not fix
this, because reads-per-molecule also depends on how many molecules each sample has -- and
lineage complexity collapses over a screen (337 UMIs/guide in the Schmierer plasmid, 122 by
day 28). Measured at matched read depth, the Schmierer high-GC deficit looked like
-0.031 (plasmid) vs -0.054 (day 4); dividing by the compression factor makes them -0.069 vs
-0.065, i.e. the same. The apparent difference was the statistic, not the chemistry.

THE FIX. Thin every sample to a common reads-per-UMI before comparing. Subsampling reads is
exactly binomial thinning of each (guide, UMI) read count, so it can be done on the table
without touching the FASTQ: draw Binomial(reads, p) per molecule and drop molecules that
reach zero -- which correctly reproduces losing them to the detection limit.

This is what "matched on molecules rather than reads" means: equalise the per-molecule
operating point, not the total.

Usage
    # thin to a target reads-per-UMI (bisects p for you)
    python3 rarefy.py $CHEM_DATA/schmierer/d4.pub.guide_umi.tsv --target-rpu 1.649

    # or a fixed fraction
    python3 rarefy.py <guide_umi.tsv> --fraction 0.5
"""
from __future__ import annotations

import argparse
import math
import random
import statistics as st
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "lib"))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

from umimodel import gc_fraction  # noqa: E402


def load(path: Path) -> dict[str, list[int]]:
    per: dict[str, list[int]] = defaultdict(list)
    with path.open() as fh:
        fh.readline()
        for line in fh:
            g, _, n = line.rstrip("\n").rsplit("\t", 2)
            per[g].append(int(n))
    return per


def thin(per: dict[str, list[int]], p: float, rng: random.Random) -> dict[str, list[int]]:
    if p >= 1.0:
        return per
    out: dict[str, list[int]] = {}
    for g, v in per.items():
        kept = []
        for n in v:
            # Binomial(n, p); exact for small n, normal approximation above 40 trials
            if n <= 40:
                k = sum(1 for _ in range(n) if rng.random() < p)
            else:
                k = int(round(rng.gauss(n * p, math.sqrt(n * p * (1 - p)))))
                k = max(0, min(n, k))
            if k:
                kept.append(k)
        if kept:
            out[g] = kept
    return out


def rpu(per: dict[str, list[int]]) -> float:
    """Median over guides of reads per observed UMI -- the operating point to match."""
    return st.median([sum(v) / len(v) for v in per.values() if v])


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("guide_umi", type=Path)
    ap.add_argument("--target-rpu", type=float, default=None,
                    help="thin until median reads-per-UMI equals this")
    ap.add_argument("--fraction", type=float, default=None)
    ap.add_argument("--min-distinct", type=int, default=20)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out-prefix", type=Path, default=None)
    a = ap.parse_args()
    rng = random.Random(a.seed)

    per = load(a.guide_umi)
    print(f"{a.guide_umi.name}: {len(per):,} guides, "
          f"{sum(len(v) for v in per.values()):,} UMIs, "
          f"{sum(sum(v) for v in per.values()):,} reads")
    print(f"  reads per UMI (median over guides): {rpu(per):.4f}")

    if a.target_rpu is not None:
        if a.target_rpu > rpu(per):
            raise SystemExit(f"target {a.target_rpu} is ABOVE this sample's "
                             f"{rpu(per):.4f}; thinning can only lower it. Pick the "
                             f"shallowest sample as the target.")
        lo, hi = 1e-4, 1.0
        for i in range(12):
            mid = 0.5 * (lo + hi)
            got = rpu(thin(per, mid, random.Random(a.seed)))
            print(f"  bisect {i + 1:2d}: p={mid:.4f} -> reads/UMI {got:.4f}")
            if got < a.target_rpu:
                lo = mid
            else:
                hi = mid
        p = 0.5 * (lo + hi)
    elif a.fraction is not None:
        p = a.fraction
    else:
        raise SystemExit("give --target-rpu or --fraction")

    out = thin(per, p, rng)
    print(f"\nkept fraction p={p:.4f}: {len(out):,} guides, "
          f"{sum(len(v) for v in out.values()):,} UMIs, "
          f"{sum(sum(v) for v in out.values()):,} reads")
    print(f"  reads per UMI now {rpu(out):.4f}")

    pref = a.out_prefix or a.guide_umi.with_name(
        a.guide_umi.name.replace(".guide_umi.tsv", "") + ".rarefied")
    pg = Path(f"{pref}.guide.tsv")
    n = 0
    with pg.open("w") as fh:
        fh.write("guide\tgc\treads\tdistinct_umi\tm1\tm2\n")
        for g, v in sorted(out.items()):
            if len(v) < a.min_distinct:
                continue
            d, r = len(v), sum(v)
            fh.write(f"{g}\t{gc_fraction(g):.4f}\t{r}\t{d}\t{r / d:.6f}\t"
                     f"{sum(x * x for x in v) / d:.6f}\n")
            n += 1
    print(f"wrote {pg}  ({n:,} guides with >= {a.min_distinct} UMIs)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
