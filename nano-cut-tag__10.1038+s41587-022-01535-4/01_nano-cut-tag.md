# nano-CUT&Tag — source notes

Defining paper: Bartosovic et al., *Nature Biotechnology* (2023), [doi:10.1038/s41587-022-01535-4](https://doi.org/10.1038/s41587-022-01535-4).

- 🟢 PMC full text, Methods and Supplementary Table 1 were read.
- 🟢 Secondary nanobody–Tn5 fusions bind mouse or rabbit primary antibodies and carry distinct MeA barcode oligos.
- 🟢 Targeted MeA tagmentation is followed by 10x Chromium linear amplification, a second standard-Tn5 tagmentation with MeB/Me-Rev, and PCR.
- 🟢 Supplementary Table 1 prints the exact representative BcdA top oligo as
  S5–`TCCACGC`–`TATAGCCT`–`GCGATCGAGGACGGC`–ME, plus MeB and
  5′-phosphorylated Tn5_Rev. Only the terminal 19-nt ME duplex is loaded on Tn5.
- 🟢 The Methods state the molecular order: deterministic MeA tagmentation, 10x
  droplet linear amplification and cell barcoding, MeB-only second tagmentation,
  then 10x scATAC PCR.
- 🟢 The published run is 36-8-48-36 (R1-I1-I2-R2), with exact custom R1 and I2
  primers and standard R2/I1.
- 🟡 The public molecule uses the BcdA branch; the source supplies BcdB/BcdC and
  staggered 0N–3N variants.
- 🔴 The complete 10x gel-bead oligo and SI-PCR primer sequences are not printed.
  Their molecular position is retained, but their bases and exact PCR landing geometry
  are not reconstructed.
