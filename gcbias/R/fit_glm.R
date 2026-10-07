#!/usr/bin/env Rscript
# Negative-binomial GLM for the GC effect. This is the preferred estimator.
#
#   reads ~ GC-band + offset(log(distinct_umi))
#
# WHY A GLM, AND WHY THIS FORM
#
# 1. It fixes the bias that broke the first analysis. Taking the median of per-guide
#    log counts is biased by -1/(2*lam*ln2), which is -0.14 at 5 UMIs per guide and
#    -0.001 at 550 -- so the bias differed between the datasets AND between GC bands
#    within a dataset, manufacturing trend. A GLM uses a log LINK on the natural count
#    scale, never a log TRANSFORM of the counts, so that class of bias does not arise.
#
# 2. log(distinct_umi) enters as a FREE COVARIATE, not an offset. An offset fixes its
#    coefficient at 1, i.e. asserts reads are exactly proportional to molecule count.
#    They are not: the fitted exponent is 1.29-1.34 on Schmierer (heavy zero-truncation
#    at ~1.6 reads/UMI makes reads grow faster than distinct UMIs) and 0.93 on
#    CRISPR-MIP. Forcing it to 1 SHRINKS the band coefficients about fourfold and
#    disagrees with the non-parametric stratified estimate; freeing it reproduces that
#    estimate. --offset re-runs the misspecified version for comparison.
#
# 3. Negative binomial, not Poisson. Fitted theta runs ~4-5 on this data, i.e. well
#    short of Poisson. Using Poisson would understate every standard error.
#
# 4. GC AS A BAND FACTOR, NOT A LINEAR TERM. On the Schmierer data a spline beats a
#    linear GC term at p = 0, and the linear slope is so unstable it changes sign
#    between replicates (-0.009, -0.052, +0.027) -- because the relationship is not
#    monotone. A band factor assumes no functional form at all and is directly
#    comparable to the binned tables in the markdown. Use --spline to check curvature.
#
# No MLE machinery beyond glm.nb is needed, and no Stan: there is no latent variable
# here once distinct_umi is used as an offset. Stan would earn its place only for a
# hierarchical model that pools across samples with the molecule count itself latent --
# worth doing after UMI error-collapse, not before.
#
# Usage:
#   Rscript fit_glm.R $CHEM_DATA/crisprmip/probe0.2_rep1.guide.tsv
#   Rscript fit_glm.R <file> --spline        # show the fitted GC curve
#   Rscript fit_glm.R <file> --offset        # misspecified version, for comparison

suppressPackageStartupMessages({library(MASS); library(data.table); library(splines)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: fit_glm.R <*.guide.tsv> [--spline] [min_distinct]")
f <- args[1]
use_spline <- "--spline" %in% args
mind <- suppressWarnings(as.integer(args[!grepl("^--", args)][2]))
if (is.na(mind)) mind <- 3L

d <- fread(f)[distinct_umi >= mind & reads > 0]
cat(sprintf("%s\n  %s guides (distinct_umi >= %d), median reads %.0f, median UMIs %.0f\n",
            basename(f), format(nrow(d), big.mark = ","), mind,
            median(d$reads), median(d$distinct_umi)))

BRK <- c(0, .45, .55, .60, .70, 1.01)
d[, band := cut(gc, BRK, right = FALSE)]
ref <- levels(d$band)[2]                      # 0.45-0.55 as reference
d[, band := relevel(factor(band), ref = ref)]

use_offset <- "--offset" %in% args
m <- if (use_offset) glm.nb(reads ~ band + offset(log(distinct_umi)), data = d) else
                     glm.nb(reads ~ band + log(distinct_umi), data = d)
if (!use_offset) {
  ex <- summary(m)$coefficients["log(distinct_umi)", ]
  cat(sprintf("\n  exponent on log(distinct_umi) = %.3f +/- %.3f  (an offset would force 1)\n",
              ex[1], ex[2]))
  if (abs(ex[1] - 1) > 5 * ex[2])
    cat("  -> significantly different from 1: an offset here would be misspecified\n")
}
co <- summary(m)$coefficients
k <- grepl("^band", rownames(co))
out <- data.table(band = sub("^band", "", rownames(co)[k]),
                  log2_est = co[k, 1] / log(2),
                  log2_se  = co[k, 2] / log(2),
                  p        = co[k, 4])
cat(sprintf("\n  NB GLM, reference band %s, theta = %.2f\n", ref, m$theta))
cat("  coefficients are log2 reads-per-molecule relative to the reference:\n\n")
print(out[, .(band, log2_est = round(log2_est, 4), log2_se = round(log2_se, 4),
              p = signif(p, 3))])

# Is a linear GC term defensible here? Report it, do not assume it.
m_lin <- glm.nb(reads ~ gc + log(distinct_umi), data = d)
m_spl <- glm.nb(reads ~ ns(gc, 4) + log(distinct_umi), data = d)
p_lrt <- anova(m_lin, m_spl)$"Pr(Chi)"[2]
cat(sprintf("\n  linear GC slope: %+.4f (log2 per unit GC %+.4f)\n",
            coef(m_lin)["gc"], coef(m_lin)["gc"] / log(2)))
cat(sprintf("  spline-vs-linear LRT p = %.3g  -> %s\n", p_lrt,
            ifelse(p_lrt < 0.01,
                   "LINEARITY REJECTED; do not quote the linear slope",
                   "linear not rejected (may simply be low power)")))

if (use_spline) {
  g <- data.table(gc = seq(0.25, 0.9, by = 0.05),
                  distinct_umi = median(d$distinct_umi))
  g[, fit := predict(m_spl, newdata = g, type = "link") / log(2)]
  g[, fit := fit - fit[which.min(abs(gc - 0.5))]]
  cat("\n  spline fit, log2 reads-per-molecule relative to GC=0.50:\n")
  print(g[, .(gc, log2 = round(fit, 4))])
}
