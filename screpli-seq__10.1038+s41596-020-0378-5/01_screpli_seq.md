# scRepli-seq — research notes

## Citation

Miura H et al. *Mapping replication timing domains genome wide in single mammalian cells with single-cell DNA replication sequencing.* Nature Protocols 15, 4058–4100 (2020). [doi:10.1038/s41596-020-0378-5](https://doi.org/10.1038/s41596-020-0378-5)

## Sources read

- 🟢 Defining Nature Protocols article, publisher PDF of record (`s41596-020-0378-5.pdf`, printed doi:10.1038/s41596-020-0378-5, Miura, Takahashi, Shibata … Hiratani & Takebayashi, Nat Protoc 15:4058–4100; hand-obtained 2026-10-17 on the institutional network, MANIFEST row). Read (`.txt` lines): Whole-genome amplification overview 570–590; Library preparation and sequencing 853–866; Comparison with Dileep & Gilbert 1160–1170; Reagents 1462–1483 (kits, TPC1/TPC2); Procedure Steps 17–79 lines 1715–2143; Step 82 trimming 2340–2425; Step 88 mapping 2500–2515; Troubleshooting 3240–3280. Steps 1–16 (cell culture, staining, micromanipulation) and the analysis Steps 89–95 were not re-read.
- 🟢 Authors' GEO deposition [GSE138634](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138634) (series PubMed 33230331; 26 GSM records). Earlier passes used it while the article was unavailable; the article now supersedes it for the method. GEO still records per-sample run lengths (80/84/127 nt, some paired-end), Ion Proton samples and SeqCap adapter use for deposited data.
- 🟢 Publisher Supplementary Information (MOESM1): Supplementary Figs 1–3 and Supplementary Note; Supplementary Software 1/2 (MOESM3/4) read-trimming and mapping procedure.
- 🟢 A later, molecularly distinct Gilbert-lab protocol (Bartlett et al., *Curr. Protoc.* 2022; [PMC8812816](https://pmc.ncbi.nlm.nih.gov/articles/PMC8812816/)) — context only, not substituted.

## 1. What it is

🟢 scRepli-seq does not label nascent DNA. It infers replication from the twofold copy-number difference between genomic domains already replicated and domains not yet replicated in one S-phase cell.

## 2. Oligos

🟢 Reagents (lines 1480–1482), verbatim with modifications:

- TPC1 primer: `5′-AATGATACGGCGACCACCGA*G-3′` (* phosphorothioate; HPLC grade; ref. 64, Quartz-seq)
- TPC2 primer: `5′-CAAGCAGAAGACGGCATACGA*G-3′` (* phosphorothioate; HPLC grade; ref. 64)

🟡 TPC1 is the first 21 nt of the Illumina P5 sequence and TPC2 the first 22 nt of P7 (scraper `contains` computation against `lib/illumina`); the page computes their landing on the assumed TruSeq-type ends.

🔴 Not printed in the article: the SeqPlex primers, the universal SEQXE adaptor sequence (Step 82 predicts one from reads), and the NEXTflex adaptor sequences/end structure.

## 3. Construct-changing steps

1. 🟢 Sort single G1 or mid-S cells (Sony SH800 single-cell mode) into 6 µl single-cell lysis and fragmentation solution (Sigma L1043 + proteinase K); 55 °C 1 h, then 99 °C 4 min, which inactivates proteinase K and fragments the gDNA (Steps 17–19; overview lines 571–573). The mid-S gate is defined in MOESM1 Supplementary Figs 1a, 3a.
2. 🟢 SeqPlex Enhanced DNA Amplification Kit (SEQXE-50XN), "based on the degenerate oligonucleotide–primed (DOP)-PCR protocol": Library Preparation Buffer 95 °C 2 min; Library Preparation Enzyme 16/24/37 °C 20 min each, 75 °C 5 min; amplification with DNA Polymerase for SeqPlex, 94 °C 2 min then 24 cycles of 94 °C 15 s / 70 °C 5 min, 70 °C 30 min (Steps 20–25).
3. 🟢 1.7× AMPure XP (Steps 26–32); NanoDrop yield distinguishes 0/1/≥2 cells (Step 33).
4. 🟢 Eco57I (Primer Removal Enzyme) digests the universal SEQXE adaptor at the ends of WGA products, 37 °C 1 h then 65 °C 20 min (Steps 34–35); the ~40-bp digested adaptors are removed with 1.8× AMPure XP (Steps 36–42), repeated if they remain (Step 43); smear 100–500 bp, peak ~150 bp (Fig. 4 legend). 🟡 Step 82 still finds a SEQXE-like sequence at some read 5′ ends, so removal is drawn as possibly incomplete.
5. 🟢 KAPA LTP (KK8232), "slightly modified … to reduce the cost" (ref. 43, Kadota 2017; the modification is not itemised): End Repair 20 °C 30 min on 20–40 ng, 1.7× AMPure (Steps 44–50); on-bead A-Tailing 30 °C 30 min, 1.8× PEG/NaCl (Steps 51–57); ligation of 5 µM NEXTflex Illumina index adaptor (NEXTflex DNA Barcodes – 48; dual-index sets allowed) with DNA Ligase 20 °C 15 min, 1.0× PEG/NaCl (Steps 58–65).
6. 🟢 Library PCR: cycle number determined by real-time PCR (KK2702, TPC primer mix, least cycles exceeding Std. 2), then KAPA HiFi HotStart ReadyMix with TPC primer mix, 98 °C 2 min, n × (98 °C 15 s, 60 °C 30 s, 72 °C 30 s), 72 °C 1 min; 5–7 cycles, published datasets 5 or 6 (Steps 66–71). 0.9× AMPure XP, repeated if ~100-bp adaptor dimers remain (Steps 72–78); library 200–700 bp (Fig. 4d).
7. 🟢 G1 cells are processed in parallel as controls (Step 18).

## 4. Final library and reads

🟡 The final product is an individually built, adaptor-indexed Illumina genomic library per cell; there is no inline cell barcode or UMI. The adaptor ends are drawn as TruSeq-type (inferred), with a 6-nt i7 window following the single-index trimming pattern in Step 82.

🟢 Sequencing: HiSeq 1500, single lane of ~40 pooled cells, "80-bp length, single read", ~4 M reads per cell (lines 860–862); Step 79 recommends >80 bp. Index read length is not stated. Step 88 maps only Read 1 (`bwa aln`/`samse`) for haplotype-resolved analysis; paired-end data may be mapped with `bwa mem`.

🟢 Step 82 prints the trimming pattern `GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG` and an example predicted SEQXE pattern `TGGTGTGTTGGGTGTGTTTCTGAAGNNNNNNNNN` (reads selected for 5′-CTGAAG-3′ in the first 50 bp). These are processing patterns, not kit oligos.

## 5. Corrections from the article (2026-10-17)

- Library PCR is performed (5–7 cycles, TPC1/TPC2), not "unstated"; the post-PCR 0.9× AMPure dimer removal is now in the graph.
- The digested SEQXE adaptors are physically removed by a 1.8× AMPure step after Eco57I; gDNA fragmentation is the 99 °C lysis step.
- The adaptor used in the Procedure is NEXTflex (single- or dual-index); SeqCap adapter kits appear only in GEO for deposited samples.
- Run layout 80-bp single read on HiSeq 1500 is from the article.

## 6. Boundary

🟢 The later DOP-PCR/barcoded pooling protocol is useful context but is molecularly distinct and is not silently substituted here.
