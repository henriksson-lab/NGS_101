# SEC-seq — source notes

Defining paper: Cheng et al., *Nature Communications* (2023), [doi:10.1038/s41467-023-39367-8](https://doi.org/10.1038/s41467-023-39367-8).

- 🟢 The PMC full text and supplementary methods were read.
- 🟢 Antibody-coated 55-µm nanovials retain cells and capture secreted IgG; PE anti-IgG followed by a barcoded anti-PE antibody reports the captured secretion in the main human B-cell experiment.
- 🟢 Cell-loaded nanovials enter a standard 10x Chromium 3′ v3.1 workflow. Transcript and antibody-derived-tag libraries are pooled 80:20 and sequenced 28-bp Read 1, 10-bp i7, 10-bp i5 and 90-bp Read 2.
- 🟢 The paper also describes a separate Raji/hybridoma compatibility experiment using TotalSeq-C streptavidin (BioLegend 405271) and anti-IgG1 (406636) with the Chromium 5′ v2 kit. Those reagents and that kit should not be conflated with the main 3′ v3.1 anti-PE readout.
- 🟡 The main experiment's commercial anti-PE barcode and complete conjugated oligo are represented by role and length because the paper does not identify their exact sequence. The 10x v3.1 bead and feature-capture oligos are also kit-defined rather than printed in this paper.
- 🟡 Raw reads are listed under GEO `GSE229042` and BioProject `PRJNA953084`. A feature-barcode consensus from those reads could help identify a used barcode after comparison with a vendor reference, but would not establish the complete antibody-conjugated oligo or its chemical modifications.
