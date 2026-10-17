# scONE-seq

**Defining paper:** Yu et al., *Science Advances* 9, eabp8901 (2023), doi:10.1126/sciadv.abp8901. Read from the PMC open-access copy PMC9812385.1 (JATS XML and PDF) and the published supplement `sciadv.abp8901_sm.pdf` (Europe PMC zip). No preprint is used.

## Sources read (round 2026-10-17-breadth-first-17)

- Article JATS `.xml.txt`: Results line 111 (workflow overview); Materials and Methods lines 211 (adaptor annealing, Tn5 loading), 215 (standard preamplification), 221 (low-volume version), 225 (Tn5 titration), 229 (library construction and sequencing).
- Supplement `.pdf.txt`: Supplementary Text "The molecular design of scONE-seq" lines 15–60; Fig. S1C legend lines 72–77; Table S2 lines 382–545; Table S3 line 732.
- Oligo scraper run over the bundle; all hits are the Table S2 rows and the two adaptor sequences repeated in Methods lines 215 and 229.

## Oligos (Table S2, 5′→3′, verbatim; spaces as printed)

🟢 One-Tn5 adaptor `[biotin]GTCTCGTGGGCTCGG TCATG NNNNNNNN AGATGTGTATAAGAGACAG`; P5-Tn5 adaptor `TCGTCGGCAGCGTC AGATGTGTATAAGAGACAG`; Mosaic `[phos]CTGTCTCTTATACACATCT`.

🟢 RT primers, mixed RT-polyT : RT-GGG : RT-TTT = 1.5 : 5 : 5 (10 µM working): `[biotin]GTCTCGTGGGCTCGG AACTT NNNNNNNN TTTTTTTTTTTTTTTTTTTTVN`, `[biotin]GTCTCGTGGGCTCGG AACTT NNNNNNNN GGG HN`, `[biotin]GTCTCGTGGGCTCGG AACTT NNNNNNNN TTT VN`.

🟢 3′ adaptor `GTCTCGTGGGCTCGG AACTT NNNNNNNN GGG HN` (RT-GGG without biotin); Amp-primer `GATGTGTGGAG GTCTCGTGGGCTCGG`.

🟢 P7 index primers `CAAGCAGAAGACGGCATACGAGAT [i7] AGATGTGTATAAGAGACAG GTCTCGTGGGCTCGG`; P5 index primers `AATGATACGGCGACCACCGAGATCTACAC [i5] TCGTCGGCAGCGTC`; P5 primer `[phos]AATGATACGGCGACCACCGA`; P7 primer `CAAGCAGAAGACGGCATACGA`.

🟢 ONE-Illu Read 1 `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`; ONE-Illu Read 2 `CCGAGCCCACGAGACCTGTCTCTTATACACATCT`; ONE-Illu I7 Index Primer `AGATGTGTATAAGAGACAGGTCTCGTGGGCTCGG`; ONE-Illu I5 Index Primer row is an asterisk ("follow the instruction provided by the manufacturers of your sequencing equipment"). ONE-MGI rows are listed for the MGI conversion route and are not modelled.

## Steps (Methods line 215, standard volume; all in one tube)

1. 🟢 Sort a cell or nucleus into 2.3 µl lysis buffer; proteinase K 55 °C 15 min, inactivate 75 °C 15 min.
2. 🟢 Tagment in 5 µl with One-Tn5 transposome (Novoprotein Tn5 loaded 24 °C 45 min with One-Tn5/Mosaic), PEG-8000, TAPS, RNase inhibitor and KAPA polymerase; 55 °C 10 min then 72 °C 10 min. 🟡 The 72 °C KAPA step is read as gap filling so that both fragment ends carry handle/DNA tag/N8/ME; the paper does not state the mechanism.
3. 🟢 Proteinase K 55 °C 15 min, 75 °C 15 min.
4. 🟢 RT with SuperScript III and the RT primer mix; 12 °C→42 °C gradient, 50 °C 50 min, 55 °C 50 min. Internal priming captures polyA and non-polyA RNA (Results line 111).
5. 🟢 Thermolabile Exo I (37 °C 20 min, 65 °C 20 min), RNase If and RNase H remove residual primers and RNA (printed "37°C for 15°C", presumably 15 min).
6. 🟢 Terminal transferase + dCTP, 37 °C 5 min; proteinase K stop.
7. 🟢 Second strand with the 3′ adaptor and KAPA: 72 °C 5 min, 10 × (48 °C 1 min, 72 °C 1 min), 72 °C 5 min. (Methods prints the 3′ adaptor as "10 mM"; Table S2 gives 100 µM stock, 10 µM working.)
8. 🟢 Exo I removes residual 3′ adaptor.
9. 🟢 One-primer co-amplification of DNA and cDNA with the Amp-primer, KAPA HotStart ReadyMix and DMSO, 18–22 cycles.
10. 🟢 AMPure XP; dilute to 0.1 ng/µl; second tagmentation with P5-Tn5 (55 °C 15 min); index PCR with P5 and the modified P7 index primers, KAPA HiFi HotStart, 10–11 cycles; pool, AMPure XP; scDASH depletion of ribosomal and mitochondrial sequences (refs 33, 58, not inspected); optional double size selection (Methods line 229).
11. 🟢 NextSeq 500 with customized sequencing primers (Table S2); Table S3 lists a High Output 150-cycle kit. Read lengths are not stated.

## Final library and reads (🟡 derived)

The Fig. S1C legend says only the 5′ (tagged) end of each preamplified fragment is kept. Assembling Table S2 oligos gives, P5→P7: P5–[i5]–`TCGTCGGCAGCGTC`–ME–insert–ME′ (DNA) or RT/3′-adaptor priming end′ (RNA)–N8′–tag′ (`CATGA` DNA / `AAGTT` RNA)–`CCGAGCCCACGAGAC`–`CTGTCTCTTATACACATCT`–[i7]–P7′. The page computes primer landing on the DNA-branch library.

## Discrepancies shown on the page

- **Read 2 / i7 labels.** The row labelled ONE-Illu Read 2 is the exact reverse complement of the row labelled ONE-Illu I7 Index Primer. On the derived library the "I7" sequence extends into the tag, UMI and insert (a Read 2 site), and the "Read 2" sequence extends into i7. The Supplementary Text says the Read 2 sequence includes the DNA and RNA barcodes, which fits the computed assignment. The page keeps the printed labels and assigns roles by computed landing. Whether the labels are a typographical swap is not established.
- **Tag length.** Results line 111 describes 6-nucleotide DNA and RNA barcodes, while Table S2 and Methods line 215 print five-base groups (`TCATG`, `AACTT`). The model uses the printed bases.

## Open questions

- i5 index-primer identity (left to the instrument manufacturer), index lengths and acquired read lengths.
- scDASH reagents and the fate of depleted molecules (cited method not inspected).
- Fate of internal P5-Tn5-only fragments from the second tagmentation (no separation described).
- MATQ-seq (ref 27), the poly(C)/3′-adaptor predecessor (ref 28) and the linked Bio-protocol procedure were not inspected.
