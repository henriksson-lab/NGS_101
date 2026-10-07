#!/usr/bin/env python3
"""Expected read counts from a lineage-UMI read-count distribution.

The question this exists to answer: for one sgRNA, is its read count higher or lower than
its *molecule* count says it should be, and does that discrepancy track guide GC%?

Naively you would divide reads by the number of distinct UMIs seen. That is wrong in two
ways that both get worse with abundance, so it manufactures exactly the trend we are
looking for:

  1. **Label collisions.** Two lineages that draw the same UMI are counted once. With a
     6-nt RSL there are only 4096 labels, and a Schmierer guide carries thousands of
     lineages, so most labels are multiply occupied.
  2. **Unseen molecules.** A molecule that happens to get zero reads is invisible, which
     biases the denominator the other way.

Both are modelled here, so that "expected reads" is a real prediction rather than a
restatement of the read count.

The model
---------
For one guide, with a label space of ``M`` possible UMIs:

* ``N`` lineages (molecules) are present; each draws a UMI uniformly, so the number of
  molecules sharing a given label is Poisson with mean ``lam = N / M``.
* Each molecule is amplified and sequenced into a read count drawn from a negative
  binomial with mean ``mu`` and shape ``r`` -- Gamma-Poisson, i.e. Poisson sequencing of a
  Gamma-distributed amplification factor. ``r -> inf`` is the Poisson (no amplification
  noise) limit; small ``r`` is a jackpotting PCR.
* A *label*'s read count is therefore the sum over however many molecules landed on it.
  Being a compound-Poisson sum of negative binomials,

      P(label gets 0 reads) = exp(-lam * (1 - p0))        with p0 = (r/(r+mu))**r

  and so, over all ``M`` labels,

      E[distinct labels observed]  D = M * (1 - exp(-lam * (1 - p0)))
      E[total reads]               R = N * mu = M * lam * mu

Those two equations are the whole model. ``D`` is what dedup gives you; ``R`` is what the
read counter gives you; the model says how they are related when nothing GC-dependent is
happening. A guide whose observed ``R`` departs from the ``R`` predicted by its ``D`` is
amplifying unusually -- and that residual is what gets regressed on GC%.

Fitting
-------
``fit_moments`` recovers ``(lam, mu, r)`` from three observables: the fraction of the label
space occupied, and the first two moments of reads per *observed* label. Writing
``q = P(0 reads) = 1 - D/M`` and taking moments of the zero-truncated distribution:

      lam * (1 - p0)              = -ln(q)                ... (1)
      lam * mu                    = m1 * (1 - q)          ... (2)
      mu * (1 + 1/r + lam)        = m2/m1 - 1             ... (3)

(3) is the second moment of the compound sum: ``E[S^2|S>0]/E[S|S>0] = 1 + mu + mu/r +
lam*mu``. The final ``lam*mu`` term is variance contributed by labels carrying more than
one molecule, and dropping it -- which is tempting, since it is absent from the
per-molecule distribution -- makes ``r`` absorb the collisions and collapse toward zero on
saturated data. Since ``lam*mu = R/M``, that term is just total reads over the label space,
so (3) stays closed-form. (2) then gives ``lam`` in terms of ``mu``, and substituting into
(1) leaves one equation in ``mu`` alone, solved by bisection. This is a
first-order (method-of-moments) fit by design: transparent, and it uses only the shape of
the UMI read-count distribution, which is the quantity we actually trust.

What actually breaks under saturation
------------------------------------
Worth being precise, because the obvious guess is wrong. The *sampling* variance of the
corrected estimator does **not** blow up as ``D -> M``; it improves, because
``lam = -ln(1 - D/M)`` grows faster than the noise in the shrinking pool of empty labels
(``molecules_rel_se`` computes this, and it falls from ~16% at 1% occupancy to ~3% at 98%).
That is a property of the estimator under a *correct* model, and it is not reassurance.

Two things do break, and both scale with read depth:

* **Naive dedup saturates.** Using ``D`` itself as the molecule count is bounded by ``M``,
  so it compresses abundant guides and nothing else. This is the trap the model exists to
  avoid, and ``simulate.py`` shows how large the induced false GC slope can be.
* **Sequencing error inflates ``D``.** An error in the UMI moves a read to a neighbouring
  label, so error reads *create* occupancy. Since ``lam`` depends on the number of empty
  labels, and that number is small when occupancy is high, a modest error rate at high depth
  drives ``D -> M`` and ``N_hat -> inf``. ``reads_before_error_fills_space`` gives the depth
  at which this takes over: for a 6-nt RSL it lands uncomfortably close to the per-guide
  depth of the Schmierer input. This -- not variance -- is the binding constraint, it is
  abundance-dependent, and abundance correlates with GC, so it confounds the exact
  comparison we are making.

The practical consequence is that guides must be filtered on occupancy, UMIs must be
error-collapsed before counting (the Schmierer authors ship a prefix-truncation script for
precisely this reason), and the surviving occupancy range has to be reported.
"""
from __future__ import annotations

import math
import random
from dataclasses import dataclass

__all__ = [
    "p0_nb", "expected_distinct", "expected_reads", "estimate_molecules",
    "molecules_rel_se", "saturation", "truncated_mean", "Fit", "fit_moments",
    "predict_reads", "log2_bias", "simulate_guide", "gc_fraction",
    "umi_error_rate", "error_inflated_distinct", "reads_before_error_fills_space",
    "chao1_labels", "good_turing_coverage", "labels_to_molecules",
    "nonparametric_molecules",
]


# --------------------------------------------------------------- the distribution itself
def p0_nb(mu: float, r: float) -> float:
    """P(a single molecule yields zero reads), negative binomial with mean mu, shape r.

    r = inf is the Poisson limit.
    """
    if mu <= 0:
        return 1.0
    if math.isinf(r):
        return math.exp(-mu)
    if r <= 0:
        return 1.0
    return math.exp(r * (math.log(r) - math.log(r + mu)))


def truncated_mean(mu: float, r: float) -> float:
    """Mean reads per molecule, conditioned on the molecule being seen at all."""
    p0 = p0_nb(mu, r)
    if p0 >= 1.0:
        return float("inf")
    return mu / (1.0 - p0)


# ------------------------------------------------------------------ forward (model -> data)
def expected_distinct(n_molecules: float, space: int, mu: float, r: float) -> float:
    """E[number of distinct UMIs observed] for n_molecules lineages in `space` labels."""
    if space <= 0:
        raise ValueError("space must be positive")
    lam = n_molecules / space
    return space * (1.0 - math.exp(-lam * (1.0 - p0_nb(mu, r))))


def expected_reads(n_molecules: float, mu: float) -> float:
    """E[total reads] -- every molecule contributes mu reads on average, seen or not."""
    return n_molecules * mu


# ------------------------------------------------------------------ inverse (data -> model)
def saturation(distinct: float, space: int) -> float:
    """Fraction of the label space occupied. Approaching 1 means dedup is meaningless."""
    return distinct / space


def estimate_molecules(distinct: float, space: int, mu: float, r: float) -> float:
    """Invert expected_distinct: how many lineages explain `distinct` observed UMIs.

    Corrects for both label collisions and molecules that drew zero reads. In the sparse,
    well-sequenced limit this collapses to `distinct`, as it must.
    """
    if distinct >= space:
        return float("inf")
    if distinct <= 0:
        return 0.0
    seen = 1.0 - p0_nb(mu, r)
    if seen <= 0:
        return float("inf")
    lam = -math.log(1.0 - distinct / space) / seen
    return lam * space


def molecules_rel_se(distinct: float, space: int, mu: float, r: float) -> float:
    """Relative standard error of estimate_molecules from sampling noise in `distinct`.

    Delta method: d(lam)/d(D) = 1 / ((1 - p0) * (space - D)), and Var(D) ~ D(space-D)/space.

    NOTE this *decreases* as D -> space, and that is not a licence to deduplicate saturated
    data. It only says the estimator is precise *if the model holds*; it says nothing about
    `distinct` being inflated by UMI sequencing error, which is the real problem at depth.
    Pair it with reads_before_error_fills_space.
    """
    if distinct <= 0 or distinct >= space:
        return float("inf")
    seen = 1.0 - p0_nb(mu, r)
    lam = -math.log(1.0 - distinct / space) / seen
    if lam <= 0:
        return float("inf")
    var_d = distinct * (space - distinct) / space
    se_lam = math.sqrt(var_d) / (seen * (space - distinct))
    return se_lam / lam


# ------------------------------------------------------------------------------ the fit
@dataclass(frozen=True)
class Fit:
    """Result of fit_moments. `ok` is False when the moments are out of model range."""
    lam: float
    mu: float
    r: float
    n_molecules: float
    ok: bool
    reason: str = ""


def fit_moments(space: int, distinct: float, m1: float, m2: float) -> Fit:
    """Recover (lam, mu, r) from the occupancy and the first two moments of reads/UMI.

    m1, m2 are E[S] and E[S^2] over *observed* labels (so every S >= 1). See the module
    docstring for the three equations; this solves them by bisection in mu.
    """
    bad = lambda why: Fit(float("nan"), float("nan"), float("nan"), float("nan"), False, why)
    if not (0 < distinct < space):
        return bad("distinct outside (0, space)")
    if m1 < 1.0:
        return bad("m1 < 1 but observed labels have >= 1 read")
    if m2 < m1 * m1:
        return bad("m2 below m1^2 -- not a valid second moment")

    q = 1.0 - distinct / space          # P(a label gets zero reads)
    neg_ln_q = -math.log(q)
    # m2/m1 - 1 = mu*(1 + 1/r) + lam*mu, and lam*mu = m1*(1-q) = reads/space, so the
    # collision contribution is observable and can be subtracted before solving for r.
    lam_mu = m1 * (1.0 - q)
    disp = m2 / m1 - 1.0 - lam_mu       # = mu * (1 + 1/r), so mu < disp strictly
    if disp <= 0:
        return bad("second moment fully explained by collisions -- r unidentifiable "
                   f"(m2/m1-1={m2 / m1 - 1:.3g}, reads/space={lam_mu:.3g})")

    def residual(mu: float) -> float:
        """Equation (1) with lam from (2) and r from (3). Zero at the solution."""
        r = mu / (disp - mu) if disp > mu else float("inf")
        lam = lam_mu / mu
        return lam * (1.0 - p0_nb(mu, r)) - neg_ln_q

    lo, hi = 1e-12 * disp, disp * (1 - 1e-12)
    f_lo, f_hi = residual(lo), residual(hi)
    if f_lo * f_hi > 0:
        # mu -> disp is the Poisson limit (r -> inf, no amplification overdispersion).
        # There the root sits exactly on the boundary rather than inside it, so a vanishing
        # residual at `hi` is a fit, not a failure.
        if abs(f_hi) <= 1e-2 * neg_ln_q:
            lam = lam_mu / hi
            return Fit(lam, hi, float("inf"), lam * space, True, "Poisson limit (r -> inf)")
        # Otherwise occupancy and shape genuinely disagree under this model -- typically
        # because `distinct` is inflated by UMI sequencing error, or the UMI pool is not
        # uniform. Both make the label space effectively smaller than `space`.
        return bad(f"no bracketed root (f={f_lo:.3g}..{f_hi:.3g}); "
                   f"occupancy {distinct / space:.3f}")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if residual(mid) * f_lo > 0:
            lo = mid
        else:
            hi = mid
    mu = 0.5 * (lo + hi)
    r = mu / (disp - mu) if disp > mu else float("inf")
    lam = lam_mu / mu
    return Fit(lam, mu, r, lam * space, True)


# ----------------------------------------------------------------- the actual comparison
def predict_reads(distinct: float, space: int, mu: float, r: float) -> float:
    """Reads this guide *should* have, given how many distinct UMIs it shows.

    `mu` and `r` are the GLOBAL amplification parameters -- fitted across all guides, so
    they carry no per-guide information. That is what makes the comparison meaningful:
    the prediction uses only the guide's molecule count, and any shortfall or excess is
    the guide amplifying differently from the library average.
    """
    return expected_reads(estimate_molecules(distinct, space, mu, r), mu)


def log2_bias(reads: float, distinct: float, space: int, mu: float, r: float) -> float:
    """log2(observed reads / predicted reads). Positive = over-amplified."""
    pred = predict_reads(distinct, space, mu, r)
    if pred <= 0 or not math.isfinite(pred) or reads <= 0:
        return float("nan")
    return math.log2(reads / pred)


# ------------------------------------- non-parametric route: no amplification model at all
# The compound-NB fit above assumes one negative binomial describes reads per molecule.
# Real PCR jackpots harder than that: on the Schmierer plasmid input the read-per-UMI
# spectrum is far more heavy-tailed than any NB that also reproduces the observed number of
# distinct UMIs, so fit_moments is rejected by the data. These estimators make no
# assumption about the amplification distribution, which is why they are the ones to
# trust when the parametric fit fails.
def chao1_labels(distinct: int, f1: int, f2: int) -> float:
    """Chao1 estimate of occupied labels, including those that got zero reads.

    f1 and f2 are the number of labels seen exactly once and exactly twice. The classic
    "unseen species" estimator: singletons tell you how much you have not seen. With no
    doubletons it falls back to the bias-corrected form.
    """
    if distinct <= 0:
        return 0.0
    if f2 > 0:
        return distinct + f1 * f1 / (2.0 * f2)
    return distinct + f1 * (f1 - 1) / 2.0


def good_turing_coverage(reads: int, f1: int) -> float:
    """Fraction of reads belonging to labels we have already seen >= twice.

    1 - f1/R. Low coverage means most of the library is represented by singletons, i.e.
    the sample is too shallow for any dedup-based statement.
    """
    if reads <= 0:
        return float("nan")
    return 1.0 - f1 / reads


def labels_to_molecules(labels: float, space: int) -> float:
    """Undo UMI collisions: occupied labels -> molecules, assuming uniform UMI draw.

    Poisson occupancy: labels = space * (1 - exp(-N/space)).
    """
    if labels <= 0:
        return 0.0
    if labels >= space:
        return float("inf")
    return -space * math.log(1.0 - labels / space)


def nonparametric_molecules(distinct: int, f1: int, f2: int, space: int) -> float:
    """Chao1 for unseen labels, then the occupancy inversion for collisions.

    This is the estimator to quote when fit_moments fails: it needs no amplification model,
    only that UMIs are drawn uniformly, and it is a lower bound rather than a point
    estimate (Chao1 under-counts when the abundance distribution is very uneven).
    """
    return labels_to_molecules(chao1_labels(distinct, f1, f2), space)


# ------------------------------------------------- sequencing error in the UMI itself
def umi_error_rate(umi_len: int, per_base_error: float) -> float:
    """P(a read's UMI carries at least one basecalling error)."""
    return 1.0 - (1.0 - per_base_error) ** umi_len


def error_inflated_distinct(true_distinct: float, reads: float, space: int,
                            umi_len: int, per_base_error: float) -> float:
    """First-order estimate of observed distinct UMIs once error reads are included.

    Error reads land on labels other than their own. Treating those as spread over the
    space, they can only add occupancy among the labels still empty, so

        D_obs ~ space - (space - D_true) * exp(-R * p_err / space)

    This is a bound-style approximation, not a correction to apply blindly: real errors land
    on *neighbouring* labels rather than uniformly, which makes the fill-in somewhat slower
    at low occupancy and faster at high. Use it to decide whether error matters at a given
    depth, and collapse UMIs properly if it does.
    """
    empty = max(0.0, space - true_distinct)
    p_err = umi_error_rate(umi_len, per_base_error)
    return space - empty * math.exp(-reads * p_err / space)


def reads_before_error_fills_space(space: int, umi_len: int,
                                   per_base_error: float) -> float:
    """Per-guide read depth at which error reads alone would cover the whole UMI space.

    Defined as R with R * p_err = space, i.e. one error read per label. Above this, the
    occupancy of a guide says more about how deeply it was sequenced than how many
    molecules it had.
    """
    p_err = umi_error_rate(umi_len, per_base_error)
    if p_err <= 0:
        return float("inf")
    return space / p_err


# --------------------------------------------------------------------------- utilities
def gc_fraction(seq: str) -> float:
    """GC content of a guide spacer, ignoring non-ACGT."""
    s = [c for c in seq.upper() if c in "ACGT"]
    if not s:
        return float("nan")
    return sum(c in "GC" for c in s) / len(s)


def simulate_guide(n_molecules: int, space: int, mu: float, r: float,
                   rng: random.Random | None = None) -> dict[str, int]:
    """Draw one guide's worth of data under the model. Returns {label: reads}, reads >= 1.

    Used by the self-test to check that the estimators recover what was put in, and by
    simulate.py to show how badly the naive ratio misbehaves under saturation.
    """
    rng = rng or random.Random(0)
    counts: dict[int, int] = {}
    for _ in range(n_molecules):
        label = rng.randrange(space)
        # NB(mean mu, shape r) = Poisson(Gamma(r, mu/r)); r = inf is plain Poisson
        rate = mu if math.isinf(r) else rng.gammavariate(r, mu / r)
        k = _poisson(rate, rng)
        if k:
            counts[label] = counts.get(label, 0) + k
    return {f"L{k}": v for k, v in counts.items()}


def _poisson(rate: float, rng: random.Random) -> int:
    """Knuth for small rate, normal approximation above 30 -- plenty for simulation."""
    if rate <= 0:
        return 0
    if rate < 30:
        el, k, p = math.exp(-rate), 0, 1.0
        while True:
            p *= rng.random()
            if p <= el:
                return k
            k += 1
    return max(0, int(round(rng.gauss(rate, math.sqrt(rate)))))
