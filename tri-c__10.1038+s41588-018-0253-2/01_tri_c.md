# Tri-C — research notes

## Citation

Oudelaar AM et al. *Single-allele chromatin interactions identify regulatory hubs in dynamic compartmentalized domains.* Nature Genetics 50, 1744–1751 (2018). [doi:10.1038/s41588-018-0253-2](https://doi.org/10.1038/s41588-018-0253-2)

## Sources read

- 🟢 The defining author manuscript and Methods in [PMC6265079](https://pmc.ncbi.nlm.nih.gov/articles/PMC6265079/).
- 🟢 The later Capture-C protocol's Tri-C module in [PMC7613269](https://pmc.ncbi.nlm.nih.gov/articles/PMC7613269/), used as corroboration.
- 🟢 The defining paper's Supplementary Tables 5 and 7 (publisher supplementary PDF), its PMC supplementary figure/table files, and the earlier bioRxiv version were checked for oligos. The cited [NG Capture-C predecessor](https://doi.org/10.1038/nmeth.3664) was checked for inherited library construction.
- 🟢 The defining article also cites the author's [2018 Tri-C Protocol Exchange method](https://doi.org/10.1038/protex.2018.113); its [2019 revision](https://doi.org/10.21203/rs.2.1650/v2) was checked independently for reagent identity and reaction order.

## 1. What it is

🟢 Tri-C targets multiway rather than only pairwise chromatin contacts. A single physical 3C concatemer must retain the viewpoint and at least two reporter fragments, so its two proximity-ligation junctions represent a concurrent single-allele contact.

## 2. Construct-changing steps

Article = defining Online Methods, "Tri-C – Experimental procedure" (PMC JATS lines 169–179). ProtEx = directly cited 2018 Protocol Exchange v1 (doi:10.1038/protex.2018.113), pp. 5–8 and rendered Tables 3 and 5 (pp. 10–11). The article says 3C libraries and capture followed its ref. 41 (NG Capture-C) "as previously described"; no step below is taken from that predecessor's own text.

1. 🟢 (Article) Fix about 10<sup>7</sup> cells in 2% formaldehyde for 10 min, quench with glycine, lyse, Dounce-homogenize, treat with 0.28% SDS then 1.67% Triton X-100.
2. 🟢 (Article) Digest with NlaIII, three 600-U additions over 16–24 h at 37 °C; heat-inactivate at 65 °C for 20 min. NlaIII was chosen because the four viewpoints lie on small NlaIII fragments (Supplementary Table 5).
3. 🟢 (Article) Proximity-ligate with 720 U T4 DNA ligase overnight at 16 °C; de-crosslink with Proteinase K overnight at 65 °C; RNase, phenol-chloroform extraction and ethanol precipitation.
4. 🟢 (Article) Covaris S220 sonication to ~450 bp to retain multiple junctions; AMPure XP clean-up and size selection at 0.7:1.
5. 🟢 (ProtEx pp. 5–6) NEBNext end repair (20 °C 30 min), dA-tailing with Klenow (37 °C 30 min), Quick T4 ligation of the kit adaptor (20 °C 15 min), then 3 µl USER enzyme added to the ligation (37 °C 15 min). The article states that the adaptors are Illumina TruSeq and that clean-ups use AMPure XP at 1.8:1. 🔴 The adaptor strands and the uracil targeted by USER are not printed in the inspected sources.
6. 🟢 (Article; ProtEx Table 3) Seven-cycle Herculase II indexing PCR with "NEB Universal primer" and "NEB Index primer" (Table 3: one pass plus "go to step 2, 6 cycles").
7. 🟢 (Article) Pool 1.5 µg of each indexed library into one multiplexed capture. (ProtEx gives 1–2 µg.)
8. 🟢 (ProtEx pp. 7–8) Round-1 capture: COT DNA plus TS-HE Universal and TS-HE Index blocking oligos; denature at 95 °C; hybridize with the biotinylated probe pool at 47 °C for 64–72 h; bind to M-270 streptavidin beads; wash, including two stringent washes at 47 °C. DNA is not eluted: 20 µl of bead suspension per reaction is amplified with KAPA HiFi HotStart ReadyMix and "POST-LM_PCR Oligo 1&2" (Table 5: 98/60/72 °C, "go to step 2, 14 cycles"). 🟡 Read in the same way as Table 3, this is 15 cycles.
9. 🟢 (ProtEx p. 8) Round-2 capture uses all amplified material (up to 2 µg), a 24-h hybridization, the same washes and the same post-capture PCR. The method states that the second capture improves enrichment 100–1000-fold.

🟡 **Published bead-ratio discrepancy.** The article prints a 0.7:1 AMPure clean-up and size selection after sonication. ProtEx step I.3 (p. 5) adds 220 µl beads to the 120-µl sonicate, about 1.8×; its note E says that all clean-ups use 1.8 volumes and that size selection is unnecessary. The graph uses the article's ratio for the published experiment, and the public page shows both values.

## 3. Final library and reads

🟢 The first erythroid experiment was sequenced on a MiSeq (v2, 250-bp paired-end reads); the second multiplexed experiment was sequenced on a NextSeq (v2, 150-bp paired-end reads). ProtEx recommends at least 300 paired-end cycles. 🔴 The index-read cycles are not stated. The molecular unit is viewpoint–reporter 1–reporter 2 (or more), not three independently observed pairwise reads.

🟢 Random sonication ends serve as unique molecular identifiers for PCR duplicate removal.

🟡 Capture-probe bases are locus-specific and therefore remain role tokens on the generic page.

🟢 Published Supplementary Table 5 supplies the genomic coordinates and 120-bp length of four Tri-C capture oligos, but does not print their base strings. Supplementary Table 7 prints two PCR primers for the separate C-Trap Nanopore validation assay; they are not Tri-C library or capture primers. The preprint numbers these as Supplementary Tables 3 and 5.

🔴 The selected 2018 TruSeq/NEBNext adapter product and sample-index pair, and the exact adapter-end ligation geometry, still need a matching historical vendor guide. The linked NG Capture-C method gives the reagent families and end-repair/dA-tailing order, but does not identify a unique Tri-C lot or index pair.

🟢 The 2018 Protocol Exchange method specifies NEBNext E6040S/L, indexed oligo sets E7500S/L and E7335S/L, end repair, dA-tailing, adapter ligation, USER treatment, then Herculase II indexing PCR. These kit identities narrow the missing vendor evidence but do not disclose the selected index oligo in each experimental library or the adapter duplex bases. The later revision repeats these reagent families and introduces Ultra II as a low-input alternative; that alternative is not evidence for the original experiment.

🟡 **Published probe-concentration discrepancy.** The 2018 Nature Genetics Methods (Tri-C experimental procedure) states 2.9 nM for the *total* equimolar pool; its directly cited 2018 Protocol Exchange method (Capture enrichment, “Biotinylated capture oligonucleotides”) specifies a working solution with *each individual oligonucleotide* at 2.9 nM. The earlier bioRxiv v1 Methods instead calls the pool 13 fmol. The latter is an amount, so it cannot be compared with concentration without the volume and number of probes. The 2019 revised method also says 2.9 nM per oligo. The model uses no concentration; the 2018 experiment's exact mixture remains uncertain.

🟡 **Published kit and dilution discrepancies.** The 2018 Protocol Exchange materials list (p. 2; extracted line 50) specifies NEBNext E6040S/L, while its library-preparation heading (p. 5; line 135) calls the procedure E7370. The 2019 revision's materials list (p. 3; line 71) and procedure heading (p. 9; lines 223–227) both specify E6040; it separately offers E7645 Ultra II for limited input (p. 9; lines 230–235). We use the 2019 revision to identify the internally consistent E6040 procedure, without assuming that every 2018 library used one kit lot. Both protocol versions also say that dissolving 10 nmol of probe in 3.46 µl yields 2.9 µM (2018 p. 7; lines 228–233; 2019 p. 14; lines 383–389). That recipe computes to about 2.9 mM, a thousandfold higher concentration. No numeric probe concentration is used in the model pending clarification.

## Bounded pass 2026-10-17 (round 5)

This was a capped pass over three originals: the article JATS Methods (lines 169–179), Protocol Exchange v1 text pp. 4–9 (extracted lines 130–325), and rendered Tables 3 and 5 (pp. 10–11). The scraper found only the C-Trap primers again.

Corrections made: added post-capture KAPA HiFi/Post-LM-PCR amplification after each capture round; separated heat inactivation, extraction, USER treatment, the AMPure clean-ups and multiplex pooling; added the hybridization blockers and conditions; named T4 ligase and Herculase II; corrected the read layout to MiSeq 250 / NextSeq 150 paired-end; added the 0.7× versus 1.8× bead discrepancy.

Deferred:
- the Protocol Exchange `supplement0.docx`;
- the NimbleGen SeqCap EZ and NEBNext historical guides, which would supply the Post-LM-PCR, blocking and adaptor bases;
- ProtEx v2 detail beyond the earlier passes;
- the Davies 2015 predecessor;
- the exact cycle-count reading;
- the index-read cycles.
