# Joung one-PCR pooled CRISPR screen readout

Joung et al. amplify an integrated sgRNA cassette directly into a sequencing-ready single-index Illumina library in one PCR. The forward primer binds hU6; the reverse primer binds the invariant sgRNA scaffold, making the geometry independent of cPPT placement.

## Source read

- 🟢 Joung et al. 2017, *Nature Protocols*, [doi:10.1038/nprot.2017.016](https://doi.org/10.1038/nprot.2017.016), especially the NGS-Lib-Fwd and NGS-Lib-Rev primer sets and the 80-cycle Read 1 plus 8-cycle index run.

🟢 Table 3 prints **ten staggered NGS-Lib-Fwd primers** and says all ten are mixed for
guide-library sequencing. It also prints **eight indexed NGS-Lib-KO-Rev primers** for
the GeCKO branch; the number used depends on sample multiplexing. The molecular page
shows one published forward/reverse member as a representative molecule, not the whole
primer mix. The SAM reverse family in the same table binds a different scaffold and is
outside this GeCKO readout.

## Screen readout steps (2026-10-17 bounded pass)

Read in the PMC author manuscript PMC5526071 (JATS text lines 856–906, 1120–1160, 2160–2214).

- 🟢 Steps 32–33: one 50 µl reaction per NGS-Lib-Fwd primer with one barcoded
  NGS-Lib-KO-Rev primer (0.25 µM each) in NEBNext High-Fidelity 2× PCR Master Mix;
  95 °C 5 min, cycles 2–13 (12 cycles) of 98 °C 20 s / 60 °C 15 s / 72 °C 15 s, 72 °C 1 min.
- 🟢 Step 58: screen gDNA uses the same reactions, scaled so all gDNA is amplified,
  up to 2.5 µg gDNA per 50 µl reaction.
- 🟢 Steps 59–64: PCRs are purified on Zymo-Spin V columns, pooled and quantified;
  Step 64 then refers to Steps 37–39. The ~260–270 bp gel extraction (Step 35) is
  written for the plasmid-library branch; the gDNA branch does not cite it.
- 🟢 Step 37: 80 cycles of Read 1 and 8 cycles of Index 1 on a MiSeq or NextSeq,
  with 5% (MiSeq) or 20% (NextSeq) PhiX.
- 🟡 The model library is computed from Fwd-1 and KO-Rev-1: {{= len(pooled_crispr.one_pcr_library().top()) }} bp
  for a non-G-starting 20-nt guide. The other forward primers change the length by their stagger.
- Open: the Discussion (line 2406) recommends titrating gDNA input and choosing the
  lowest visible-band cycle number, so the 12 cycles are the printed default, not a
  fixed screen value. Supplementary code and sample-specific index assignment are not checked.

## Distinguishing feature

One PCR simultaneously enriches the integrated guide and adds the complete P5/P7 library arms. It is not interchangeable with workflows where PCR1 first enriches the locus and PCR2 later adds indexed adapters.
