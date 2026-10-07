#!/usr/bin/env Rscript
# Qualitative plots -- look at these before believing any slope.
#
# Usage:
#   Rscript plot_qualitative.R $CHEM_DATA/schmierer/input.25000000.guide.tsv
#
# Writes PDFs next to the input. Nothing is written into the repo.

suppressPackageStartupMessages({library(data.table); library(ggplot2)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: plot_qualitative.R <*.guide.tsv> [umi_len]")
f <- args[1]
umi_len <- if (length(args) > 1) as.integer(args[2]) else 6L
space <- 4^umi_len

g <- fread(f)
g[, occupancy := distinct_umi / space]
g[, reads_per_umi := reads / distinct_umi]
stem <- sub("\\.guide\\.tsv$", "", f)

theme_set(theme_bw(base_size = 10))

# 1 -- is the UMI space saturated? This gates everything else.
p1 <- ggplot(g, aes(occupancy)) +
  geom_histogram(bins = 60, fill = "grey30") +
  geom_vline(xintercept = 0.9, colour = "firebrick", linetype = 2) +
  labs(title = sprintf("UMI-space occupancy per guide (space = %s)",
                       format(space, big.mark = ",")),
       subtitle = "right of the dashed line, distinct-UMI counts cannot be trusted",
       x = "distinct UMIs / label space", y = "guides")

# 2 -- the confound: does GC track abundance? If yes, naive ratios are unsafe.
p2 <- ggplot(g, aes(gc, reads)) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), colour = "steelblue") +
  scale_y_log10() +
  labs(title = "Abundance vs GC% -- the confound",
       subtitle = "a slope here means a reads/UMI trend may be biology, not amplification",
       x = "guide GC fraction", y = "reads (log10)")

# 3 -- the naive statistic, which saturation distorts
p3 <- ggplot(g, aes(gc, reads_per_umi)) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), colour = "firebrick") +
  labs(title = "NAIVE reads per distinct UMI vs GC%",
       subtitle = "compare against the model statistic in plot_gcbias.R before interpreting",
       x = "guide GC fraction", y = "reads / distinct UMI")

# 4 -- occupancy against abundance: shows the saturation mechanism directly
p4 <- ggplot(g, aes(reads, occupancy)) +
  geom_point(alpha = 0.15, size = 0.5) +
  geom_hline(yintercept = 1, colour = "firebrick", linetype = 2) +
  scale_x_log10() +
  labs(title = "Why the naive ratio fails: occupancy is a function of depth",
       subtitle = "the ceiling at 1.0 compresses exactly the abundant guides",
       x = "reads (log10)", y = "occupancy")

pdf(sprintf("%s.qualitative.pdf", stem), width = 7, height = 5)
for (p in list(p1, p2, p3, p4)) print(p)
dev.off()

# 5 -- the per-UMI read spectrum, if extract.py produced it
spec_f <- sprintf("%s.umi_spectrum.tsv", stem)
if (file.exists(spec_f)) {
  s <- fread(spec_f)
  p5 <- ggplot(s[reads_per_umi <= 40], aes(reads_per_umi, n_umis)) +
    geom_col(fill = "grey30") + scale_y_log10() +
    labs(title = "Per-UMI read-count spectrum (the shape the model fits)",
         subtitle = "a long jackpot tail means small r: amplification is overdispersed",
         x = "reads per UMI", y = "number of UMIs (log10)")
  pdf(sprintf("%s.umi_spectrum.pdf", stem), width = 7, height = 4); print(p5); dev.off()
  cat(sprintf("wrote %s.umi_spectrum.pdf\n", stem))
}

# GC quintile summary, if present
q_f <- sprintf("%s.gc_quintiles.tsv", stem)
if (file.exists(q_f)) print(fread(q_f))

cat(sprintf("wrote %s.qualitative.pdf\n", stem))
