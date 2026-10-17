# MUSIC — source notes

Defining paper: Zhao et al., *Nature* (2024), [doi:10.1038/s41586-024-07239-w](https://doi.org/10.1038/s41586-024-07239-w).

- 🟢 PMC full text, Methods and supplementary notes were read.
- Round 2026-10-17-breadth-first-13 (bounded pass): published JATS
  `s41586-024-07239-w_PMC11023937.xml.txt` Critical reagents lines 169–203 and the
  procedure lines 221–258 (RNA-linker ligation → sequencing); Supplementary Note 2
  (MOESM1 text lines 100–200); Supplementary Table 2 sheet (`export.xlsx.txt`,
  Linkers and Adaptors rows 1–15 and member 1 of each barcode sheet). The 10x
  CG000315 Rev E guide remains a manual download (curl DNS error) and was not read.
- 🟢 RNA linker: `5′-OH-CGAGGAGCGCTTNNNNNrArUrArGrCrArUrUrGrC-OH-3′`; the five random DNA bases are a UMI.
- 🟢 The Y-shaped DNA linker has printed top and bottom strands and a 5-nt UMI. A single-base sticky end ligates to dA-tailed genomic DNA.

🟢 **Published linker-length discrepancy.** The defining paper's Critical reagents DNA-linker paragraph (published PDF text lines 1811–1819) says the hybridized top and bottom strands make a **15-bp** duplex and a single 3′-dT overhang. Its later DNA-linker ligation paragraph (lines 1991–1995) calls those same printed strands a **14-bp** duplex. Supplementary Table 2, Linkers and Adaptors rows 3–4, gives the strands verbatim; their complementary overlap is 15 bp with the stated one-base overhang. The model uses those published sequences and their 15-bp pairing. Why the ligation paragraph counts 14 bp remains unresolved.
- 🟢 Three sequential 14-bp barcode sets, each with 96 members, encode the cell. The complex identity combines a 16-bp 10x GEM barcode with an 8-bp i7; the adjacent 12-bp 10x UMI deduplicates molecules.
- 🟢 The run is 28-bp Read 1, 8-bp index and 150-bp Read 2. Read 2 contains cell barcodes, analyte linker and RNA or DNA insert.
- 🟢 The authors' public oligo spreadsheet was downloaded as XLSX. It prints the complete universal PCR adapter, 24 indexed PCR adapters, both linker strands and all three 96-member barcode sets. The page shows the universal primer, index A1 and member 1 of each barcode round as copyable representatives.
- 🟢 The same ordering workbook specifies a terminal phosphorothioate bond
  (`...C*T`) on the universal and indexed PCR adapters, `/5OH/` on the
  DNA-linker bottom and barcode-3 bottom strands, and `/5Phos/` on the other
  displayed barcode arms. Those modifications are retained in the page's
  copyable oligo list; the plain DNA constants remain suitable for base-pairing.
- 🟢 **Published reaction order** (Methods lines 221–243): T4 PNK 5′ phosphorylation of
  RNA → RNA linker with T4 RNA ligase 1 (15% PEG 8000, 25 °C 2 h then 16 °C overnight) →
  HpyCH4V digestion (37 °C 3 h) → dA-tailing with Klenow fragment → DNA linker with
  Instant Sticky-end Ligase (20 °C 6 h) → T4 PNK phosphorylation of the linker 5′ ends →
  three 96-well split-pool barcode ligations (overnight 20 °C, EDTA/EGTA quench) → T4 PNK
  3′ dephosphorylation → E. coli poly(A) polymerase tailing → 10-µm strainer, 5,000
  nuclei → Covaris sonication → 10x Chromium Next GEM 3′ (PN-1000268) GEM RT → eight
  aliquots, proteinase K reverse crosslinking, Monarch RNA columns → Bst 3.0 55 °C 15 min
  → 1.8× RNAClean → NEBNext Ultra II Q5 PCR with Universal Adaptor and one index adaptor
  per aliquot, 13–14 cycles → 1.2× AMPure, pool → 4% E-gel 300 bp–1.2 kb, Monarch gel
  extraction. The earlier model ligated the DNA linker before the RNA linker and lacked
  the digestion, tailing, phosphatase, sonication, crosslink-reversal and cleanup steps;
  this is a repository modelling correction.
- 🟡 RNA-linker junction (derived from Methods lines 171 and 223): the linker is
  5′-OH DNA / 3′-OH RNA and the nuclear RNA 5′ ends are phosphorylated first, so T4 RNA
  ligase 1 joins the linker 3′ end to the RNA 5′ end; the linker's own 5′ end becomes a
  donor only after the later PNK step.
- 🟢 Supplementary Note 2: RNA cDNA = cell barcode + RNA linker + RNA insert + 10X
  barcode; DNA = cell barcode + DNA-linker(bt) + DNA insert + DNA-linker(top) + 10X
  barcode. The page now draws separate RNA- and DNA-derived library molecules; linker
  regions are role placeholders sized from the printed oligos.
- 🟢 Repository transcription correction: the DNA linker top strand has **30** A after NB
  (Methods line 177 and Supplementary Table 2 row DNA_Linker_Top_PolA; consistent with the
  stated 36-nt top ssDNA). The model previously had 31.
- Open: the purpose of the Bst 3.0 incubation is not stated in the inspected Methods;
  the 10x partial Read 1 bases (CG000315) are unread, so PCR primer sites are not drawn
  on a pre-PCR template; barcode overhang junctions are not drawn.
- 🟢 The final construct order is derived from those primer sites: Read 1 reports the 16-base GEM barcode and 12-base UMI, whereas Read 2 enters barcode round 3, then rounds 2 and 1, the analyte linker and finally the RNA or DNA insert.
