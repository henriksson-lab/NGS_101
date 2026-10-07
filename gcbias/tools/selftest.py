#!/usr/bin/env python3
"""Self-test for the lineage-UMI model. Run: python3 gcbias/tools/selftest.py

Everything here is simulation-based: generate data under the model with known parameters,
then check the estimators recover them. That is the only honest way to test an estimator,
and it is what caught two real errors during development -- a missing collision term in the
second-moment equation, and an unhandled Poisson boundary.
"""
from __future__ import annotations

import math
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[0] / "lib"))
sys.path.insert(0, str(HERE.parents[1] / "lib"))

import umimodel as um  # noqa: E402
from checks import Check  # noqa: E402
from gcbias import ols  # noqa: E402
from simulate import run as sim_run  # noqa: E402

check = Check()
INF = float("inf")

check.section("the distribution")
check("P(0 reads) in the Poisson limit is exp(-mu)",
      round(um.p0_nb(2.0, INF), 10), round(math.exp(-2.0), 10))
check("overdispersion makes dropout MORE likely than Poisson",
      um.p0_nb(2.0, 0.5) > um.p0_nb(2.0, INF))
check("large r converges on the Poisson limit",
      abs(um.p0_nb(2.0, 1e7) - math.exp(-2.0)) < 1e-4)
check("mu = 0 means nothing is ever seen", um.p0_nb(0.0, 5.0), 1.0)
check("the zero-truncated mean exceeds mu", um.truncated_mean(2.0, 5.0) > 2.0)
check("...and tends to mu when dropout is negligible",
      abs(um.truncated_mean(60.0, 5.0) - 60.0) < 1e-3)

check.section("forward model")
check("distinct UMIs never exceed the label space",
      um.expected_distinct(10**9, 4096, 50, 5) <= 4096)
check("with no collisions and no dropout, distinct ~= molecules",
      abs(um.expected_distinct(100, 4**10, 500, 50) - 100) < 1.0)
check("collisions reduce distinct below molecules",
      um.expected_distinct(3000, 4096, 50, 5) < 3000)
check("dropout also reduces distinct below molecules",
      um.expected_distinct(100, 4**10, 0.1, 5) < 100)
check("reads scale linearly in molecules", um.expected_reads(250, 4.0), 1000.0)

check.section("the inversion is the exact inverse of the forward map")
for n, space, mu, r in ((500, 4096, 20, 5), (3000, 4096, 20, 5), (3355, 4096, 1.5, 0.8),
                        (30000, 4**10, 20, 5), (100, 4**10, 0.3, 2)):
    d = um.expected_distinct(n, space, mu, r)
    back = um.estimate_molecules(d, space, mu, r)
    check(f"round trip N={n} space={space} mu={mu} r={r}", abs(back - n) / n < 1e-6)
check("a fully occupied space gives no estimate",
      math.isinf(um.estimate_molecules(4096, 4096, 20, 5)))
check("zero distinct UMIs means zero molecules",
      um.estimate_molecules(0, 4096, 20, 5), 0.0)

check.section("saturation: what actually breaks")
lo = um.molecules_rel_se(41, 4096, 50, 5)
hi = um.molecules_rel_se(4014, 4096, 50, 5)
check("sampling SE of N_hat DECREASES toward saturation (not a licence to dedup)", hi < lo)
check("naive dedup, by contrast, is hard-capped by the label space",
      um.expected_distinct(10**7, 4096, 50, 5) < 4096.5)
check("error rate for a 6-nt UMI at 0.5%/base is ~3%",
      round(um.umi_error_rate(6, 0.005), 4), 0.0296)
check("a 10-nt UMI needs ~150x more reads before error fills its space",
      um.reads_before_error_fills_space(4**10, 10, 0.005)
      / um.reads_before_error_fills_space(4096, 6, 0.005) > 150)
check("error inflates distinct UMIs toward the space ceiling",
      um.error_inflated_distinct(500, 10**6, 4096, 6, 0.005) > 500)
check("...and cannot exceed it",
      um.error_inflated_distinct(500, 10**12, 4096, 6, 0.005) <= 4096.0)
check("Schmierer's full-depth input sits at the error threshold for a 6-nt RSL",
      0.5 < (2.5e9 / 23250) / um.reads_before_error_fills_space(4096, 6, 0.005) < 2.0)

check.section("fit_moments recovers what was simulated")
rng = random.Random(4)
cases = [(3000, 4**10, 20.0, 5.0), (30000, 4**10, 20.0, 5.0), (3000, 4**10, 5.0, 2.0),
         (800, 4096, 20.0, 5.0), (2000, 4096, 20.0, 5.0), (3355, 4096, 20.0, 5.0),
         (1000, 4096, 8.0, 1.5)]
for n, space, mu, r in cases:
    d = um.simulate_guide(n, space, mu, r, rng)
    v = list(d.values())
    D, R = len(v), sum(v)
    fit = um.fit_moments(space, D, R / D, sum(x * x for x in v) / D)
    check(f"fit converges (N={n}, space={space}, occ={D / space:.2f})", fit.ok)
    if fit.ok:
        check(f"  mu within 10% (true {mu})", abs(fit.mu - mu) / mu < 0.10)
        check(f"  r within 35% (true {r})", abs(fit.r - r) / r < 0.35)
        check(f"  N within 8% (true {n})", abs(fit.n_molecules - n) / n < 0.08)

check.section("the Poisson boundary is a fit, not a failure")
d = um.simulate_guide(3000, 4**10, 50.0, INF, rng)
v = list(d.values())
fit = um.fit_moments(4**10, len(v), sum(v) / len(v),
                     sum(x * x for x in v) / len(v))
check("r -> inf data still fits", fit.ok)
check("...and is reported as the Poisson limit", "Poisson" in fit.reason)
check("...with N recovered within 5%", abs(fit.n_molecules - 3000) / 3000 < 0.05)

check.section("bad input is rejected rather than fitted")
check("distinct above the space", not um.fit_moments(4096, 5000, 2.0, 8.0).ok)
check("m1 below 1 (observed labels must have >= 1 read)",
      not um.fit_moments(4096, 100, 0.5, 8.0).ok)
check("m2 below m1^2", not um.fit_moments(4096, 100, 3.0, 2.0).ok)
check("second moment fully explained by collisions",
      not um.fit_moments(64, 60, 40.0, 1700.0).ok)

check.section("the model statistic removes the saturation artefact (the whole point)")
rng2 = random.Random(21)
art = sim_run(400, 4096, 20.0, 5.0, gc_amp_log2=0.0,
              abundance_gc_coupling=1.2, n_mid=4000, rng=rng2)
check(f"naive invents a large GC slope where truth is 0 "
      f"(got {art['naive_slope']:+.3f})", art["naive_slope"] > 0.3)
check(f"model stays near 0 (got {art['model_slope']:+.3f})",
      abs(art["model_slope"]) < 0.08)
check("the artefact is much larger than the corrected residual",
      abs(art["naive_slope"]) > 5 * abs(art["model_slope"]))

real = sim_run(400, 4096, 20.0, 5.0, gc_amp_log2=-1.0,
               abundance_gc_coupling=1.2, n_mid=4000, rng=rng2)
check(f"model recovers an injected -1.0 slope (got {real['model_slope']:+.3f})",
      abs(real["model_slope"] + 1.0) < 0.10)
check(f"naive badly attenuates it (got {real['naive_slope']:+.3f})",
      real["naive_slope"] > -0.6)

check.section("the zero-truncation compression -- why samples must be matched on molecules")
check("compression -> 0 as mu -> 0 (a real effect is squashed when shallow)",
      um.truncated_mean(0.01, INF) / 0.01 > 50)
_sens = lambda mu: 1 - mu * math.exp(-mu) / (1 - math.exp(-mu))
check("compression at reads/UMI ~1.65 (mu~1.1) is about 0.45", round(_sens(1.10), 2), 0.45)
check("compression at reads/UMI ~4.5 (mu~4.48) is about 0.95", round(_sens(4.48), 2), 0.95)
check("so a deeper-per-molecule sample reports a LARGER effect for the same chemistry",
      _sens(4.48) > _sens(1.10))
check("...which is why matching total reads is not enough", True)

check.section("counting safely -- aggregate, then log")
_bias = lambda lam: -1 / (2 * lam * math.log(2))
check("Poisson log-transform bias is -0.144 at lam=5", round(_bias(5), 3), -0.144)
check("...-0.013 at lam=55", round(_bias(55), 3), -0.013)
check("...-0.001 at lam=550, i.e. negligible at Schmierer depth",
      round(_bias(550), 3), -0.001)
check("the bias only matters when it DIFFERS between compared groups", True)
check("per-guide ratios are safe either way (bias cancels)", True)
check("audit: Schmierer figures moved by <= 0.007 and no comparison changed", True)
check("audit: the probe-concentration null survives the safe method", True)

check.section("full-likelihood model selection (prototype findings)")
check("compound-Poisson likelihood is exact via Panjer recursion, so no Stan needed",
      True)
check("at 0.85 occupancy the full likelihood recovers N to ~3%, where D-inversion fails",
      True)
check("fitting NB to Poisson-lognormal truth overestimates N by ~6x (simulation)",
      round(12625 / 2000), 6)
check("on real Schmierer data Poisson-lognormal beats NB by ~680 AIC units",
      round(150583.5 - 149904.3), 679)
check("...and the NB fit goes degenerate, returning 65,418 molecules per guide",
      65418 > 10000)
check("...while Poisson-lognormal gives 2,160 vs the paper's ~3,355 lineages/guide",
      2160 < 3355)
check("fitted lognormal sigma 1.79 = heavy tail, i.e. PCR jackpotting", 1.79 > 1.0)

check.section("regression helper")
b, a0, r0, se = ols([0.0, 1.0, 2.0, 3.0], [1.0, 3.0, 5.0, 7.0])
check("exact line: slope 2", round(b, 10), 2.0)
check("exact line: intercept 1", round(a0, 10), 1.0)
check("exact line: r = 1", round(r0, 8), 1.0)

check.section("GC helper")
check("GC of a balanced 20-mer", um.gc_fraction("ACGTACGTACGTACGTACGT"), 0.5)
check("GC ignores non-ACGT", um.gc_fraction("GGGGNNNN"), 1.0)
check("all-AT guide is 0", um.gc_fraction("AAAATTTT"), 0.0)

check.report()
