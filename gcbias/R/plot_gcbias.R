#!/usr/bin/env Rscript
# The comparison that answers the reviewers: model-corrected bias vs GC%, with the naive
# statistic overlaid so the difference between them is visible rather than asserted.
#
# Usage:
#   Rscript plot_gcbias.R $CHEM_DATA/schmierer/input.25000000.gcbias.tsv

suppressPackageStartupMessages({library(data.table); library(ggplot2); library(mgcv)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: plot_gcbias.R <*.gcbias.tsv>")
f <- args[1]
d <- fread(f)
stem <- sub("\\.gcbias\\.tsv$", "", f)
theme_set(theme_bw(base_size = 10))

long <- rbind(
  data.table(gc = d$gc, value = d$bias_log2,  statistic = "model: log2(obs / predicted)"),
  data.table(gc = d$gc, value = d$naive_log2, statistic = "naive: log2(reads / distinct UMI)")
)

p1 <- ggplot(long, aes(gc, value, colour = statistic, fill = statistic)) +
  geom_point(alpha = 0.12, size = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs")) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_colour_manual(values = c("steelblue", "firebrick")) +
  scale_fill_manual(values = c("steelblue", "firebrick")) +
  labs(title = "GC bias: model-corrected vs naive",
       subtitle = "divergence between the curves is saturation artefact, not chemistry",
       x = "guide GC fraction", y = "log2 deviation (centred)") +
  theme(legend.position = "bottom", legend.direction = "vertical")

# Binned medians -- robust, and the form to quote qualitatively.
d[, gc_bin := cut(gc, breaks = seq(0, 1, by = 0.05), include.lowest = TRUE)]
b <- d[, .(n = .N,
           model = median(bias_log2), naive = median(naive_log2),
           lo = quantile(bias_log2, 0.25), hi = quantile(bias_log2, 0.75),
           gc = mean(gc)), by = gc_bin][n >= 10][order(gc)]

p2 <- ggplot(b, aes(gc, model)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.2, fill = "steelblue") +
  geom_line(colour = "steelblue") + geom_point(colour = "steelblue") +
  geom_line(aes(y = naive), colour = "firebrick", linetype = 2) +
  geom_hline(yintercept = 0, linetype = 3) +
  labs(title = "Binned median bias (blue = model, IQR ribbon; red dashed = naive)",
       subtitle = "look for movement in the end bins, not just the overall slope",
       x = "guide GC fraction", y = "median log2 deviation")

# Occupancy check: the model statistic must not depend on occupancy. If it does, the
# correction has not fully removed the saturation and the GC slope is still suspect.
p3 <- ggplot(d, aes(occupancy, bias_log2)) +
  geom_point(alpha = 0.12, size = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), colour = "steelblue") +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(title = "Residual dependence on UMI occupancy (should be flat)",
       subtitle = "a trend here invalidates the GC slope -- tighten --max-occupancy",
       x = "occupancy", y = "model log2 deviation")

pdf(sprintf("%s.gcbias.pdf", stem), width = 7, height = 5)
for (p in list(p1, p2, p3)) print(p)
dev.off()
cat(sprintf("wrote %s.gcbias.pdf\n", stem))

cat("\nbinned medians:\n"); print(b)
cat("\nlinear slopes (log2 per unit GC fraction):\n")
for (col in c("bias_log2", "naive_log2")) {
  m <- lm(as.formula(sprintf("%s ~ gc", col)), data = d)
  s <- summary(m)$coefficients["gc", ]
  cat(sprintf("  %-12s slope %+0.4f  se %0.4f  p %.3g\n", col, s[1], s[2], s[4]))
}
cat("\nis the GC effect non-linear (do the extremes differ from the bulk)?\n")
gm <- gam(bias_log2 ~ s(gc, bs = "cs"), data = d)
print(summary(gm)$s.table)
