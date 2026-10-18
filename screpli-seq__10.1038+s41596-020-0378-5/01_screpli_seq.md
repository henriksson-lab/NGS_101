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

🟡 TPC1 is the first 21 nt of the Illumina P5 sequence and TPC2 the first 22 nt of P7 (scraper `contains` computation against `lib/illumina`); this identifies sequence overlap only. Native NEXTflex landing sites and first-cycle product geometry remain unassigned; the graph does not use a completed canonical amplicon as its PCR substrate.

The exact synthesized SeqPlex primer species, universal SEQXE adaptor strands and native NEXTflex end structure remain unassigned in the bounded scope. Step 82 predicts a SEQXE-like pattern from reads; it does not establish synthesized reagent identity. This is not a publication-wide absence claim.

## 3. Construct-changing steps

1. 🟢 Sort single G1 or mid-S cells (Sony SH800 single-cell mode) into 6 µl single-cell lysis and fragmentation solution (Sigma L1043 + proteinase K); 55 °C 1 h, then 99 °C 4 min, which inactivates proteinase K and fragments the gDNA (Steps 17–19; overview lines 571–573). The mid-S gate is defined in MOESM1 Supplementary Figs 1a, 3a.
2. 🟢 SeqPlex Enhanced DNA Amplification Kit (SEQXE-50XN), "based on the degenerate oligonucleotide–primed (DOP)-PCR protocol": Library Preparation Buffer 95 °C 2 min; Library Preparation Enzyme 16/24/37 °C 20 min each, 75 °C 5 min; amplification with DNA Polymerase for SeqPlex, 94 °C 2 min then 24 cycles of 94 °C 15 s / 70 °C 5 min, 70 °C 30 min (Steps 20–25).
3. 🟢 1.7× AMPure XP (Steps 26–32); NanoDrop yield distinguishes 0/1/≥2 cells (Step 33).
4. 🟢 Eco57I (Primer Removal Enzyme) digests the universal SEQXE adaptor at the ends of WGA products, 37 °C 1 h then 65 °C 20 min (Steps 34–35); the ~40-bp digested adaptors are removed with 1.8× AMPure XP (Steps 36–42), repeated if they remain (Step 43); smear 100–500 bp, peak ~150 bp (Fig. 4 legend). 🟡 Step 82 still finds a SEQXE-like sequence at some read 5′ ends, so removal is drawn as possibly incomplete.
5. 🟢 KAPA LTP (KK8232), "slightly modified … to reduce the cost" (ref. 43, Kadota 2017; the modification is not itemised): End Repair 20 °C 30 min on 20–40 ng, 1.7× AMPure (Steps 44–50); on-bead A-Tailing 30 °C 30 min, 1.8× PEG/NaCl (Steps 51–57); ligation of 5 µM NEXTflex Illumina index adaptor (NEXTflex DNA Barcodes – 48; dual-index sets allowed) with DNA Ligase 20 °C 15 min, 1.0× PEG/NaCl (Steps 58–65).
6. 🟢 Library PCR: cycle number determined by real-time PCR (KK2702, TPC primer mix, least cycles exceeding Std. 2), then KAPA HiFi HotStart ReadyMix with TPC primer mix, 98 °C 2 min, n × (98 °C 15 s, 60 °C 30 s, 72 °C 30 s), 72 °C 1 min; 5–7 cycles, published datasets 5 or 6 (Steps 66–71). 0.9× AMPure XP, repeated if ~100-bp adaptor dimers remain (Steps 72–78); library 200–700 bp (Fig. 4d).
7. 🟢 G1 cells are processed in parallel as controls (Step 18).

## 4. Final library and reads

🟢 The source constructs an individual adaptor-indexed library per cell. The graph selects the NEXTflex DNA Barcodes – 48 single-index variant named in Step 58; the source also allows separate dual-index and unique-dual variants. Native adapter strands, index width/member and primer landing sites remain unassigned. No inline cell barcode, UMI, canonical TruSeq shell or fixed six-base index is added. The single-index trimming example does not establish an acquired index-cycle length.

🟢 Sequencing: HiSeq 1500, single lane of ~40 pooled cells, "80-bp length, single read", ~4 M reads per cell (lines 860–862); Step 79 recommends >80 bp. Index read length is not stated. Step 88 maps only Read 1 (`bwa aln`/`samse`) for haplotype-resolved analysis; paired-end data may be mapped with `bwa mem`.

🟢 Step 82 prints the trimming pattern `GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG` and an example predicted SEQXE pattern `TGGTGTGTTGGGTGTGTTTCTGAAGNNNNNNNNN` (reads selected for 5′-CTGAAG-3′ in the first 50 bp). These are processing patterns, not kit oligos.

## 5. Corrections from the article (2026-10-17)

- Library PCR is performed (5–7 cycles, TPC1/TPC2), not "unstated"; the post-PCR 0.9× AMPure dimer removal is now in the graph.
- The digested SEQXE adaptors are physically removed by a 1.8× AMPure step after Eco57I; gDNA fragmentation is the 99 °C lysis step.
- The adaptor used in the Procedure is NEXTflex (single- or dual-index); SeqCap adapter kits appear only in GEO for deposited samples.
- Run layout 80-bp single read on HiSeq 1500 is from the article.

## 6. Boundary

🟢 The later DOP-PCR/barcoded pooling protocol is useful context but is molecularly distinct and is not silently substituted here.

## Bounded native-graph migration (2026-10-18, batch 12)

One original inspected: published article `s41596-020-0378-5.pdf.txt`, WGA overview 570–590, Library preparation and sequencing 853–866, Reagents 1462–1483, Procedure Steps 17–79 (1715–2143), and adaptor trimming/prediction 2340–2385. The mandatory bundle scraper returned ten hits in four files; relevant TPC1/TPC2 bases and terminal phosphorothioate bonds were verified at article lines 1480–1483. The false equipment-model hit at line 1533 is not an oligo. Later-protocol and software hits were classified by identity and not adopted as SeqPlex reagent bases; the article supplies the active trimming-pattern evidence in this pass.

Removed shared conformation/replication recipe imports and the inferred canonical final-library shell. Every S-phase and G1 lineage now declares its own native states and physical allocations. WGA PCR receives the actual pre-amplification population, while diagnostic and bulk library PCR receive distinct native ligated aliquots. Original genomic DNA, WGA daughters and later TPC-initiated DNA are independent populations. Physical TPC primers retain their terminal phosphorothioate bond after extension; copying a modified template does not duplicate its backbone modification.

Source allocations are explicit: 500–800 ng of purified WGA DNA enters Eco57I digestion (Step 34); 20–40 ng enters end repair (Step 44); 1.5 µl of the 11-µl ligated eluate enters diagnostic real-time PCR and 7 µl enters bulk PCR (Steps 65–71). Remaining material has an unresolved subsequent disposition rather than being reused in both reactions. Diagnostic PCR is terminal QC, not a library precursor returned to the bulk reaction. Unspecified NanoDrop/MultiNA/quantification aliquots are not invented.

The graph retains each positive source cleanup: 1.7× AMPure after WGA, 1.8× after Eco57I, 1.7× after end repair, on-bead A-tailing with 1.8× PEG/NaCl, 1.0× PEG/NaCl after ligation and 0.9× AMPure after library PCR. Released SEQXE fragments and genomic polymers are separate digestion outputs; short-fragment depletion is physical, while later SEQXE trimming and G1 normalization are computational. Residual templates, primers or adaptors within a retained fraction have unresolved component-specific recovery. Source repeat-cleanup criteria remain attached to the operations; no particular sample is asserted to have required repeats.

Deferred: applicable SeqPlex generation/manual and native ordered primer species; NEXTflex DNA Barcodes – 48 native physical strands/end geometry, cut coordinates, complete first-cycle PCR products, index member/width and index cycles; component-specific cleanup recoveries; analytical QC sample volumes; complete membership/volumes of the ~40-cell sequencing pools; linked Kadota 2017/Sasagawa 2013 originals and other article sections not reread here. Native kit priming is explicitly unresolved rather than promoted to a verified transition. The existing missing-information audit status remains unchanged.
