# Drop-BS — source notes

Defining paper: Fu et al., *Nature Communications* (2023), [doi:10.1038/s41467-023-40411-w](https://doi.org/10.1038/s41467-023-40411-w).

- 🟢 PMC full text and supplementary protocol diagrams were read.
- 🟢 MNase fragments DNA in single-nucleus droplets. Dielectrophoretic fusion pairs these with bead droplets containing end-repair/ligation reagents.
- 🟢 UV releases bead oligos and ligation covalently adds a 15-nt cell barcode before bisulfite conversion.
- 🟢 Supplementary Figure 1 prints the ligation-end arm
  `GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG` and the complementary oligo as
  `3′-/3SpC3/TCTACACATATTCTCTGTC/5Phos/-5′`; the latter is stored 5′→3′ in the
  model and annealed base-by-base before the ligation is drawn.
- 🟢 Barcoded DNA is bisulfite-converted in droplets, random-primed with S5–ME–N9, and completed by i5/i7 indexed PCR.
- 🟢 Supplementary Figure 1 also prints the intervening pre-PCR oligos. The P5 transposon primer is `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`; the P7 transposon BS primer is `GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGTAAGTAGAAGATGGTATATGAGAT` when its diagrammed 3′→5′ strand is rewritten 5′→3′. These were omitted from the earlier page's ordering list. The N5XX and Ad2.XX indexed PCR primers are diagrammed as families with unspecified i5/i7 index members, so the selected exact index oligos remain open.
- 🟢 The final diagram in the supplement supports a Nextera-compatible sequencing endpoint.
