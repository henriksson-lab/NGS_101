# Drop-BS — source notes

Defining paper: Zhang et al., *Nature Communications* (2023), [doi:10.1038/s41467-023-40411-w](https://doi.org/10.1038/s41467-023-40411-w).

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

## Round 2026-10-17-breadth-first-09 bounded pass

Read: published PMC10400590 JATS Methods lines 189–219 (bead design and complement annealing, droplet lysis/MNase, barcoding by fusion, droplet bisulfite conversion, random priming, pre-amplification, indexing PCR, quantification and sequencing); Springer MOESM1 Supplementary Fig. 1 text lines 29–130. The bioRxiv preprint was not used (published version first).

- 🟢 Bead oligo (Fig. S1): 5′-[photocleavable linker]-`CAAGCAGAAGACGGCATACGAGAT`-D15-`GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG`-3′ (ChemGenes lot L142464-1). The 19-nt complement (IDT) is annealed at 85 °C for 5 min with slow cooling.
- 🟢 Droplet lysis buffer contains MNase and CaCl2: ice 10 min, RT 15 min, then 65 °C 30 min. Fusion uses 300 Vp-p at 200 kHz. The bead/enzyme droplet holds Fast-Link ligase and End-It; the ligation mix holds EGTA, dNTPs and ATP. Then 365 nm UV for 8 min, 16 °C overnight, 70 °C 15 min, perfluorooctanol, proteinase K 65 °C 2 h, and SPRIselect at 1.4× then 1×.
- 🟢 Bisulfite: 95 °C denaturation, then MethylEdge reagent in droplets at 95 °C 5 min and 54 °C 1 h. Desulphonation/column cleanup follows the kit, then 1.1× SPRI.
- 🟢 Random priming uses Klenow exo− (Enzymatics), slow ramps to 37 °C, then 60 min, followed by 1.1× SPRI, denaturation and 1.1× SPRI. Pre-PCR uses the P5 transposon and P7 transposon BS primers with KAPA HiFi HotStart Uracil+ for 5 cycles. Indexing PCR uses N5XX/Ad2.XX with EvaGreen monitoring, then 0.85× SPRI. Sequencing is NovaSeq 6000 S4, 150 bp paired-end.
- 🟡 Derived in the model: the synthetic bead oligo is unmethylated, so every bead C reads as T after conversion. The A/G/T barcode is unchanged. The converted ligation-arm complement ({{= drop_bs.CONVERTED_ARM_RC }}) and converted PCR-site complement ({{= drop_bs.CONVERTED_PCR_SITE_RC }}) match the copy strands printed in Fig. S1. The P7 transposon BS primer equals the Nextera Read 2 arm plus the converted bead-site complement, and construction fails otherwise. Read 2 therefore reads 24 nt of converted bead site, the 15-nt barcode and 34 nt of converted arm before genomic sequence.
- 🟡 Published-symbol discrepancy, shown on the page: Methods and Fig. S1 define the barcode as A/G/T (D15), but Fig. S1 labels the copied strand V15. The complement of D is H (C/A/T). The model uses H; no printed base is affected.
- 🔴 Selected i5/i7 index members and index lengths are not printed in the sources read. The page draws both index windows without a length, and the read table reports the index reads as of unstated length. The ChemGenes lot specification was not available.

## Batch02 native-graph migration (2026-10-18)

The local endpoint is now assembled from the actual source N5XX/Ad2.XX primer families, the source P5/P7 pre-PCR primers, and the random-copy product. Index widths and selected members remain unassigned; canonical Nextera sequence constants supply only the source-identified components. The actual pre-PCR template is single-strand S5–ME–N9/random-copy DNA. The second PCR receives the preceding five-cycle pre-PCR product, including its newly installed Read2 arm, rather than the expected final indexed library. N9 is a random-priming coordinate, not a UMI. R2 traverses the converted24nt bead site before the15nt cell barcode, matching the source’s R2positions25–39 extraction. Biological PE150 acquisition is retained; index acquisition is explicitly unassigned in the read table.

The19nt blocked complement pairs flush with the bead oligo3′ terminus: a blunt duplex tail, not a ligation splint spanning bead and genome. The source /5Phos/ and /3SpC3/ modifiers are preserved. Failed droplet fusion is no longer labelled physically discarded, and original converted templates persist as an unresolved class after random extension. Physical SPRI/column/oil separation remains distinct from chemical transformation.

This bounded pass read two originals: published PMC10400590 JATS Methods191–229 (beads, droplet fusion, recovery, bisulfite, random priming, both PCRs and PE150/barcode analysis), and published MOESM1 FigS1 sequence/strand panels. Mandatory scraper pointers were verified there; preprint sources and large derived methylation tables were not used. The published D15-versus-V15 copy-symbol discrepancy remains visible, with H15 used for the actual complement of the A/G/T barcode alphabet. Detailed first-cycle copying, original-template recovery, selected index members/widths and exact index acquisition remain deferred. This is not a full review or absence acceptance.
