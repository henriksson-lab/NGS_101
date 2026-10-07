#!/usr/bin/env Rscript
# Expected vs actual read count per sgRNA, coloured by GC%.
#
# Expected = (molecules estimated from the lineage-UMI distribution) x (one global
# reads-per-molecule). So the identity line is the reference: a guide sitting off it is
# amplifying differently from the library average, which is the GC question.
#
# Usage:
#   Rscript plot_expected_vs_actual.R $CHEM_DATA/schmierer/input.collapsed.gcbias.tsv
#   Rscript plot_expected_vs_actual.R <file> nb      # use the NB fit instead of Chao1
#
# A WORD ON WHAT IS VISIBLE. Per-guide scatter is ~0.30 log2 while the whole GC effect is
# ~0.05 log2, so the raw scatter (panel 1) cannot show the effect by eye no matter how it is
# coloured -- roughly a 6:1 noise-to-signal ratio. Panels 2-4 bin it, which is the only way
# the colour becomes informative. Panel 1 is included because it is the honest picture of
# how good the prediction is per guide.

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: plot_expected_vs_actual.R <*.gcbias.tsv> [nb|chao1]")
f <- args[1]
which_pred <- if (length(args) > 1) args[2] else "chao1"

d <- fread(f)
stem <- sub("\\.gcbias\\.tsv$", "", f)

if (which_pred == "nb" || !"np_pred_reads" %in% names(d)) {
  d[, expected := pred_reads]
  pred_lab <- "NB compound model"
} else {
  d[, expected := np_pred_reads]
  pred_lab <- "Chao1 molecules x global reads/molecule"
}
d <- d[is.finite(expected) & expected > 0 & reads > 0]
d[, ratio := reads / expected]
d[, lr := log2(ratio)]

# Discrete GC bands, chosen on quantiles so each holds a comparable number of guides.
brks <- unique(quantile(d$gc, probs = seq(0, 1, length.out = 6), na.rm = TRUE))
d[, gc_band := cut(gc, breaks = brks, include.lowest = TRUE, dig.lab = 2)]

theme_set(theme_bw(base_size = 10))
# Colour limits from the data's own spread (2nd-98th percentile), not a round guess:
# most guides sit in 0.45-0.70, and fixed wide limits crush them all into one hue.
gclim <- as.numeric(quantile(d$gc, c(0.02, 0.98), na.rm = TRUE))
gcs <- scale_colour_viridis_c(option = "plasma", name = "guide GC", limits = gclim,
                              oob = squish)

rng <- range(c(d$expected, d$reads))

# --- 1. the plot as asked: expected vs actual, coloured by GC ----------------------------
p1 <- ggplot(d, aes(expected, reads, colour = gc)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey20", linewidth = 0.4) +
  geom_point(alpha = 0.35, size = 0.45) +
  gcs +
  scale_x_log10(limits = rng) + scale_y_log10(limits = rng) +
  coord_fixed() +
  labs(title = "Expected vs actual read count per sgRNA",
       subtitle = sprintf("expected = %s;  grey line = identity;  n = %s guides",
                          pred_lab, format(nrow(d), big.mark = ",")),
       x = "expected reads (from lineage-UMI molecule count)", y = "actual reads")

# --- 2. same data, GC as bands with a fit per band ---------------------------------------
# Makes the GC dependence legible: parallel lines offset from identity = a GC effect.
p2 <- ggplot(d, aes(expected, reads, colour = gc_band)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey20", linewidth = 0.4) +
  geom_point(alpha = 0.12, size = 0.35) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, formula = y ~ x) +
  scale_colour_viridis_d(option = "plasma", end = 0.9, name = "GC band") +
  scale_x_log10(limits = rng) + scale_y_log10(limits = rng) +
  coord_fixed() +
  labs(title = "Expected vs actual, split into GC bands",
       subtitle = "a band sitting below identity yields fewer reads than its molecules predict",
       x = "expected reads", y = "actual reads")

# --- 3. the ratio, which is where the effect is actually visible -------------------------
p3 <- ggplot(d, aes(expected, ratio, colour = gc)) +
  geom_hline(yintercept = 1, colour = "grey20", linewidth = 0.4) +
  geom_point(alpha = 0.3, size = 0.45) +
  gcs +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Actual / expected vs expected",
       subtitle = "flat in x means the prediction is unbiased with respect to abundance",
       x = "expected reads", y = "actual / expected")

# --- 4. binned: median actual vs median expected within GC x abundance cells -------------
# This is the panel to read. Binning collapses the 6:1 noise so the colour separates.
d[, ab_bin := cut(expected, breaks = quantile(expected, seq(0, 1, length.out = 9)),
                  include.lowest = TRUE, labels = FALSE)]
agg <- d[, .(n = .N, exp_med = as.numeric(median(expected)),
             act_med = as.numeric(median(reads)), lr_med = as.numeric(median(lr)),
             lr_se = sd(lr) / sqrt(.N), gc_med = median(gc)),
         by = .(gc_band, ab_bin)][n >= 30]

p4 <- ggplot(agg, aes(exp_med, act_med, colour = gc_band)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey20", linewidth = 0.4) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.6) +
  scale_colour_viridis_d(option = "plasma", end = 0.9, name = "GC band") +
  scale_x_log10() + scale_y_log10() + coord_fixed() +
  labs(title = "Binned medians: expected vs actual, by GC band",
       subtitle = "each point is >=30 guides; separation between lines is the GC effect",
       x = "median expected reads", y = "median actual reads")

# --- 5. the quantitative readout: offset from identity per GC band -----------------------
band <- d[, .(n = .N, lr = as.numeric(median(lr)),
              lo = as.numeric(median(lr) - 1.96 * sd(lr) / sqrt(.N)),
              hi = as.numeric(median(lr) + 1.96 * sd(lr) / sqrt(.N)),
              gc = as.numeric(median(gc))), by = gc_band][order(gc)]

p5 <- ggplot(band, aes(gc, lr)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey30") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.008, colour = "grey40") +
  geom_point(aes(colour = gc), size = 3) + gcs +
  labs(title = "Deviation from expectation by GC band",
       subtitle = "log2(actual/expected), median +/- 95% CI of the median",
       x = "median guide GC of band", y = "log2(actual / expected)")

pdf(sprintf("%s.expected_vs_actual.pdf", stem), width = 7, height = 6)
for (p in list(p1, p2, p3, p4, p5)) print(p)
dev.off()
cat(sprintf("wrote %s.expected_vs_actual.pdf  (5 panels)\n\n", stem))

cat("per-GC-band deviation from expectation:\n")
print(band[, .(gc_band, n, gc = round(gc, 3), log2_ratio = round(lr, 4),
               ci_lo = round(lo, 4), ci_hi = round(hi, 4))])

cat("\nnoise check -- why panel 1 cannot show this by eye:\n")
cat(sprintf("  per-guide sd of log2(actual/expected) : %.3f\n", sd(d$lr)))
cat(sprintf("  total GC effect across all bands      : %.3f\n",
            diff(range(band$lr))))
cat(sprintf("  noise : signal                        : %.1f : 1\n",
            sd(d$lr) / diff(range(band$lr))))

cat("\ncorrelation of actual with expected:\n")
cat(sprintf("  Pearson  (log10) %.4f\n", cor(log10(d$expected), log10(d$reads))))
cat(sprintf("  Spearman         %.4f\n",
            cor(d$expected, d$reads, method = "spearman")))

# --- abundance control -------------------------------------------------------------------
# GC correlates with depth, and the ratio is not perfectly flat in depth, so the GC bands
# could partly be an abundance effect. Stratify and check the GC pattern survives inside
# each abundance quartile.
cat("\nabundance confound check:\n")
cat(sprintf("  Spearman(log2 ratio, expected) = %+.4f  (want ~0)\n",
            cor(d$expected, d$lr, method = "spearman")))
cat(sprintf("  Spearman(GC, expected)         = %+.4f\n",
            cor(d$gc, d$expected, method = "spearman")))

d[, ab_q := cut(expected, quantile(expected, seq(0, 1, length.out = 5)),
                include.lowest = TRUE, labels = c("Q1 low", "Q2", "Q3", "Q4 high"))]
cat("\nlog2(actual/expected) by GC band x abundance quartile:\n")
print(dcast(d[, .(v = round(as.numeric(median(lr)), 4)), by = .(gc_band, ab_q)],
            gc_band ~ ab_q, value.var = "v"))

cat("\nGC effect relative to the second (mid-GC) band, within each abundance quartile:\n")
for (q in levels(d$ab_q)) {
  s2 <- d[ab_q == q, .(v = as.numeric(median(lr))), by = gc_band][order(gc_band)]
  cat(sprintf("  %-8s lowGC-mid = %+.4f   highGC-mid = %+.4f\n",
              q, s2$v[1] - s2$v[2], s2$v[nrow(s2)] - s2$v[2]))
}

p6 <- ggplot(d[, .(lr = as.numeric(median(lr)), n = .N), by = .(gc_band, ab_q)][n >= 30],
             aes(gc_band, lr, colour = ab_q, group = ab_q)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey30") +
  geom_line() + geom_point(size = 2) +
  scale_colour_viridis_d(option = "viridis", end = 0.85, name = "abundance") +
  labs(title = "GC effect within abundance quartiles (confound control)",
       subtitle = "the GC pattern should repeat inside each quartile if it is real",
       x = "GC band", y = "median log2(actual / expected)") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

pdf(sprintf("%s.expected_vs_actual_control.pdf", stem), width = 7, height = 4.5)
print(p6); dev.off()
cat(sprintf("\nwrote %s.expected_vs_actual_control.pdf\n", stem))

# PNG previews of the two panels most worth eyeballing
png(sprintf("%s.preview_scatter.png", stem), width = 1500, height = 1400, res = 180)
print(p1); dev.off()
png(sprintf("%s.preview_binned.png", stem), width = 1500, height = 1400, res = 180)
print(p4); dev.off()
png(sprintf("%s.preview_band.png", stem), width = 1500, height = 1000, res = 180)
print(p5); dev.off()
