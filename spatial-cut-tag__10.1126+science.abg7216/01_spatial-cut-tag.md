# Spatial-CUT&Tag / hsrChST-seq — source notes

## What it is

Spatial-CUT&Tag is the published name of the method introduced as
**high-spatial-resolution chromatin modification state profiling by sequencing
(hsrChST-seq)** in the authors' preprint. It is one protocol, not a family: antibody-tethered
pA–Tn5 tags histone-mark-associated DNA in a fixed tissue section, then two perpendicular
microfluidic ligations add an A×B spatial address.

## Sources read

- 🟢 Deng et al., *Science* (2022),
  [doi:10.1126/science.abg7216](https://doi.org/10.1126/science.abg7216), is the defining
  published paper. 🟢 Its Supplementary Materials (author-manuscript DOCX
  `EMS145971-supplement-Supplementary_Material.docx`, PMC7612972; hand-downloaded by the
  maintainer 2026-10-17; title and author list match) were read: Materials and Methods
  "DNA barcodes and other key reagents", "Transposome preparation", "Spatial-CUT&Tag" and
  "Data preprocessing" (DOCX text lines 902–932) and Tables S1–S3 (lines 1375–1528). This
  is now the authority for the page's oligos and conditions.
- 🟢 Deng et al., bioRxiv (2021),
  [doi:10.1101/2021.03.11.434985](https://doi.org/10.1101/2021.03.11.434985), supplies the
  original hsrChST-seq name and describes antibody binding, pA–Tn5 activation, A and B
  ligations, crosslink reversal and library PCR.
- 🟢 Farzad et al., *Nature Protocols* (2024),
  [doi:10.1038/s41596-024-01013-y](https://doi.org/10.1038/s41596-024-01013-y), is the
  authors' detailed protocol for Spatial-CUT&Tag and Spatial-ATAC-seq. The MOESM2
  workbook sheets are unlabelled by assay: Supp_Table 1 prints Tn5MErev, a Tn5ME-A,
  Tn5ME-B, both ligation linkers, the custom Read 1 primer and N501/N701–N705; Supp_Table 2
  prints the transposome/linker oligos without a Tn5ME-A plus New-P5-SMART and Ad2.1–2.8;
  Supp_Table 3 prints the 50 + 50 spatial barcodes. Read in the 2026-10-17 pass: JATS
  Reagents (lines 225–260), Reagent setup (410–444), Steps 16–24 (516–560), CUT&Tag
  Steps 36–60 (600–660) and library PCR/cleanup Steps 105–121 (850–1000).
- 🟢 The authors' `dyxmvp/spatial-CUT-Tag` processing repository was read to verify the raw
  read coordinates: barcode B is raw Read 1 bases 1–8, barcode A is bases 39–46, and raw
  Read 2 is genomic sequence.

- 🟢 GEO series [GSE165217](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165217),
  the Science paper's data deposit, read as its NCBI GEO FTP SOFT family file in the
  2026-10-17 follow-up (12 GSM records with identical protocol fields). Its extraction
  protocol summarizes the chemistry (light formaldehyde fixation; primary then secondary
  antibody; Mg2+-activated pA–Tn5 inserting "adapters containing a ligation linker"; A and B
  in-situ ligations; crosslink reversal and PCR) and states HiSeq "pair-end 150 bp mode with
  custom read 1 primer" on an Illumina HiSeq 4000. Its data-processing field says Read 1
  was first filtered by "two constant linker sequences (linker 1 and linker 2)" before
  barcodes A and B and the genomic sequence were separated. No oligo sequences are printed.

The 2021 bioRxiv PDF appends Supplementary Fig. S1 and other figure captions, but no
ordered oligo table or full wet-lab recipe. Its Data and materials availability section
links `dyxmvp/hsrChST-seq`, which did not resolve when checked; the accessible
`dyxmvp/spatial-CUT-Tag` repository may be a successor, but equivalence is unverified.
The 2024 supplementary PDF MOESM1 contains only a gene-score mapping figure; the
ordered oligos below are in the separate MOESM2 workbook. The 2022 Science
Supplementary Materials (now read) print the same Tn5MErev, Tn5ME-B, ligation linkers,
New-P5-SMART, Ad2.1–Ad2.8, custom Read 1 primer and all 50 + 50 barcodes as the 2024
tables (compared base for base); only Tn5ME-A differs (below).

## Molecular reconstruction

🟢 Deng et al. 2022 Supplementary Materials (lines 909–929), as now drawn: unloaded pA–Tn5
(Diagenode C01070002) loaded per the manufacturer with Tn5ME-A
`5′-/5Phos/CATCGGCGTACGACTAGATGTGTATAAGAGACAG-3′`, Tn5ME-B and Tn5MErev; section fixed
with 0.2% formaldehyde 5 min, 1.25 M glycine; NP40–digitonin permeabilization 5 min;
primary antibody 1:50 overnight at 4 °C; secondary antibody 1:50, 30 min RT, 5-min
wash-buffer wash; pA–Tn5 complex 1:100 in 300-wash buffer 1 h RT, 5-min 300-wash; 10 mM
MgCl2 in 300-wash buffer 37 °C 1 h; 40 mM EDTA 5 min; NEBuffer 3.1 wash. Barcodes A
(A1–A50) annealed to linker 1 and barcodes B (B1–B50) annealed to linker 2 (2X annealing
buffer); per channel 2 µL ligation mix (T4 DNA ligase, T4 ligase buffer, Triton X-100) +
2 µL NEBuffer 3.1 + 1 µL 25 µM barcode, 37 °C 30 min each; NEBuffer 3.1 or DPBS washes.
Lysis in 0.1% SDS/10 mM TAPS at 60 °C 2 h, then 65 °C 1 h with rotation (no proteinase K
named). PCR per 5 µL lysate: 15 µL 0.67% Triton X-100, 2 µL 10 µM new P5 primer, 2 µL
10 µM i7 primer, 25 µL NEBNext High-Fidelity 2X master mix; 58 °C 5 min, 72 °C 5 min,
98 °C 30 s, 12 × (98 °C 10 s, 60 °C 10 s), 72 °C 1 min; 1.3× AMPure XP; Bioanalyzer HS;
HiSeq 4000 PE150 with custom Read 1 primer.

🟢 (historical, later implementation) Farzad et al. 2024 CUT&Tag branch: NP40–digitonin permeabilization; primary antibody
(H3K27me3, H3K4me3 or H3K27ac, 1:50, 2 h RT or 4 °C overnight); secondary antibody (1:50,
30 min) and three 150-wash washes; pA–Tn5 diluted 1:100 for 1 h at RT, then two 300-wash
washes; activation in tagmentation buffer at 37 °C for 1 h; 40 mM EDTA stop; barcode A
and barcode B in-tissue ligations; lysis; PCR (58 °C 5 min, 72 °C 5 min, 98 °C 30 s, then
12 cycles) and a 1× (65 µl) AMPure XP cleanup. The procedure names its PCR primers only
as "universal or barcoded i5/i7 primer"; 🟡 the model uses New-P5-SMART and Ad2, the only
printed primers that land on the barcode-B and S7 ends. The 2022 methods name a "new P5
PCR primer" and "i7 primers" and Table S1 prints New-P5-SMART and Ad2.1–Ad2.8. Where the
2022 and 2024 conditions differ (primary-antibody time, wash counts, 1.3× vs 1× AMPure),
the page now follows 2022.

🟡 The splint geometry is computed from the printed oligos (the module raises if it does
not hold). Ligation linker 2 is the reverse complement of the barcode-B 3′ end followed
by the 5′-phosphorylated barcode-A 5′ end. Ligation linker 1 is the reverse complement of
the barcode-A 3′ mosaic end followed by `CATCGGCGTACGACT`: the adapter that pA–Tn5
transfers onto genomic DNA must therefore begin with that 15-nt extension, followed by a
mosaic end, for barcode A to be ligated. 🟢 The 2022 Tn5ME-A supplies exactly that
5′-phosphorylated extension (the module raises if it does not), so the segments are no
longer drawn as inferred.

🟢 Published discrepancy (on the page): Farzad 2024 Supp_Table 1 prints Tn5ME-A as
`/5Phos/TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`, whose extension has no complementary
stretch in linker 1. The model uses the 2022 sequence from the defining publication;
whether the 2024 row is a table error or reflects other loading oligos is not stated.
🟡 Which Tn5ME pairing each transposome carries is not stated; the drawn A–B fragment is
the end class that both ligations and the PCR require.

🟡 Model correction (2026-10-17): the page previously used the shared
`spatial_epigenome` construct, which placed linker 1's reverse complement directly
against genomic DNA and omitted the transferred adapter strand. The module now builds
its construct locally.

🟢 The 50 first-barcode sequences and 50 second-barcode sequences are public. The page
shows one representative oligo from each plate and records the two 8-base regions as
separate parts of one combinatorial spatial coordinate.

## Read layout

🟢 The authors' repository and the 2021 GEO deposit state HiSeq 4000 paired-end 150
sequencing; GEO adds that a custom Read 1 primer was used for the deposited (2022-era)
data. Raw Read 1 carries barcodes A and B and raw Read 2 genomic sequence. GEO names
the primer but not its sequence. 🟢 The 2022 Supplementary Materials (line 929) state the same run and print the
custom Read 1 primer in Table S1; "Data preprocessing" (line 932) filters Read 1 on linker
1 and linker 2 and converts it to Cell Ranger ATAC format. 🟡 (computed) The custom Read 1
primer enters barcode B (cycles 1–8), linker 2 and barcode A (cycles 39–46), then the
barcode-A mosaic end, the Tn5ME-A extension and mosaic end, so genomic DNA starts at
cycle {{= spatial_cut_tag._GENOMIC_START }}. The Nextera-compatible index primer reads the
i7 sample index, and the Nextera Read 2 primer enters genomic DNA from the S7 end.

## Open questions

- Meaning of the 2024 Supp_Table 1 Tn5ME-A row; Diagenode C01070002 loading guidelines
  (ratio of Tn5ME-A and Tn5ME-B duplexes) not inspected.
- Which Ad2 index each sample used, and the acquired index cycles.

## Follow-up pass 2026-10-17 (round deferred-followup-02)

Routes tried for the Science supplementary methods (EMS145971 DOCX): PMC bin URL (PMC
"Preparing to download" proof-of-work page, not solved by script); legacy
ncbi.nlm.nih.gov/pmc bin URL (404); PMC Cloud bucket listing for PMC7612972 (XML, TXT and
JSON only, licence TDM, no supplement); Europe PMC REST supplementaryFiles ("not open
access"), Europe PMC bin and fulltextRepo URLs (403 challenge / 500); science.org article
and supplement URLs (Cloudflare challenge); Unpaywall/OpenAlex (only the bioRxiv v1
preprint, which carries figure captions but no Materials and Methods, and the PMC author
manuscript); bioRxiv API (only v1). Manual URL remains
<https://pmc.ncbi.nlm.nih.gov/articles/instance/7612972/bin/EMS145971-supplement-Supplementary_Material.docx>.
The GEO deposit GSE165217 was read instead (above); it confirms the custom Read 1 primer
and instrument but supplies no oligos.

## Follow-up pass 2026-10-17 (round deferred-followup-03)

The defining Science Supplementary Materials DOCX is now on disk (hand-placed; MANIFEST
row records the PMC bin URL; title "Spatial-CUT&Tag: Spatially Resolved Chromatin
Modification Profiling at Cellular Level", Deng … Fan — identity matches). Read: DOCX text
lines 887–932 (Materials and Methods through data preprocessing) and 1375–1528 (Tables
S1–S3). The `.txt` twin carries embedded EndNote base64; the inspected passages are plain
text. Scraper: 24 hits; Tn5ME-A (line 912) found only here. Corrections: Tn5ME-A from the
2022 supplement replaces the inferred adapter extension; the discrepancy is now 2022 vs
2024 Tn5ME-A; conditions follow 2022 (overnight primary, single washes, SDS/TAPS lysis,
NEBNext HF master mix, 1.3× AMPure); whitelists cite 2022 Tables S1–S2; Read 1 genomic
start computed. Selftest gains a 2022 Tn5ME-A transcription check (skips if absent).
Deferred: Diagenode C01070002 loading guidelines, the 2024 Tn5ME-A explanation, Farzad
NIHMS bin files, Ad2 per-sample assignment and acquired index cycles. No
missing-information claim is accepted.
