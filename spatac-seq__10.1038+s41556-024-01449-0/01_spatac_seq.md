# SPATAC-seq — source notes

## What it is

SPATAC-seq is the **split-pool ligation-based assay for transposase-accessible chromatin
using sequencing**. Fixed nuclei receive a two-ended first barcode during custom-Tn5
transposition, then one barcode is ligated to each end in two successive 48-well
pool-and-split rounds. Indexed PCR creates sublibraries and completes the Illumina arms.

This method is unrelated to Spatial-ATAC-seq despite the similar spelling.

## Sources read

- 🟢 Sun et al., *Nature Cell Biology* (2024),
  [doi:10.1038/s41556-024-01449-0](https://doi.org/10.1038/s41556-024-01449-0), is the
  defining paper. Figure 1 gives the molecular architecture and the four indexing rounds.
  The full article PDF with Methods (Nature Cell Biology 26:1187–1199) was placed by hand
  through institutional access on 2026-10-17 and its Methods are now read (first reader;
  see "Published Methods" below). An independent second reader (2026-10-17) re-read the
  Methods and the rendered Extended Data Fig. 1a; see "Second reading" below.
- 🟢 Supplementary Table 1 prints every Tn5, ligation-splint, round-2, round-3, i5 and i7
  oligo and identifies the `5′P` and `3′ddC` modifications.
- 🟢 The authors' [Lan-lab/SPATAC-seq](https://github.com/Lan-lab/SPATAC-seq) repository
  supplies the workflow drawing, the 110,592-entry 32-nt valid-barcode list and executable
  read preprocessing. Its script extracts round-2 at Read 1 cycles 1–8, the forward half
  of round-1 at Read 1 cycles 39–46, round-3 at Read 2 cycles 1–8 and the reverse half of
  round-1 at Read 2 cycles 24–31. Genomic sequence begins at Read 1 cycle 66 and Read 2
  cycle 51.
- 🟢 GEO series [GSE215901](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE215901)
  (species mixing) and [GSE243256](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243256)
  (zebrafish atlas), read as the GEO FTP SOFT family files in the 2026-10-17 follow-up.
  GSE215901 states that 48 indexed Tn5 transposomes were assembled and that the
  species-mixing run used only the first eight, giving eight sublibraries; its library
  protocol field only says "see description of our paper". Neither record describes the
  ligation rounds or PCR. (Corrected 2026-10-17: an earlier version of this note said the
  GEO record describes the two ligation rounds and indexed PCR; it was written while the
  record was unread behind reCAPTCHA.)

## Molecular reconstruction

🟢 `Tn5_F_D1` and `Tn5_R_D1` each carry an 8-nt round-1 barcode between a 15-nt leader
and the 19-nt mosaic end. The two halves are retained as distinct structured parts of the
same combinatorial cell barcode. The authors use 48 specified forward/reverse pairs rather
than all 64 possible combinations of the eight oligos on each side.

🟢 The round-2 oligo is
`ACGACGCTCTTCCGATCT-[8 nt]-ATCCACGTGCTTGAG`. Its 3′ end ligates to the
5′-phosphorylated forward round-1 adapter. `Tn5_F_LA` is a 23-nt splint: its reverse
complement is exactly the 15-nt round-2 bridge followed by the first eight bases of the
round-1 leader. The page therefore computes and marks that ligation boundary.

🟢 The round-3 oligo is `S7-mosaic end-[8 nt]`. Its 3′ end ligates to the
5′-phosphorylated reverse round-1 adapter. Each `Tn5_R_LA-*` splint contains eight bases
complementary to the round-1 leader, the reverse complement of the corresponding
round-3 barcode and seven bases of the round-3 mosaic end. This second ligation boundary
is likewise anchored to named segments.

🟢 Round-4 primers add P5, i5 and a complete TruSeq Read 1 site on the forward end, and
P7, i7 and S7 on the reverse end. The result intentionally combines a TruSeq Read 1 entry
with Nextera Index 1 and Read 2 entries.

## Read layout

🟢 The authors' processing recipe uses raw Read 1 bases 66–115 and raw Read 2 bases
51–100 as 50-nt genomic mates. The discarded prefixes are molecular information, not
arbitrary trimming: Read 1 is `round2(8) + bridge(30) + round1-forward(8) + ME(19)`, and
Read 2 is `round3(8) + leader(15) + round1-reverse(8) + ME(19)`.

🟢 Both GEO records' data-processing fields describe the genomic mates as raw Read 1
positions 66–121 and raw Read 2 positions 51–101. The start positions agree with the code;
the end positions (121/101) differ from the code's 115/100 spans. 🟡 These are trimming
choices, not chemistry; the page uses the code's spans and shows both statements.

🟢 Acquired cycles: the Reporting Summary (MOESM2) and the GSE243256 series overall
design both state paired-end 150-bp reads, so the graph and page carry 150/8/8/150
acquired cycles (2026-10-17 shared-code follow-up; previously the page showed the code's
115/100 spans as cycle counts). The code's 115/100 spans are kept as
`PROCESSING_SPANS` in `tools/chem.py` and the selftest checks the barcode coordinates
against them. The construct, rather than total cycle counts, is the invariant: longer
reads continue farther into genomic DNA without changing the barcode positions. 🟡 On
the page's short representative 44-nt insert a 150-cycle read runs through into the
opposite barcode block; that read-through is a property of the placeholder insert
length, not of typical fragments.

## Published Methods (first reading, 2026-10-17)

Read in `s41556-024-01449-0.pdf.txt` lines 3323–3559 (Methods: nuclei isolation and
fixation, SPATAC-seq method, transposome assembly, ligation-oligo preparation,
transposition, hybridization and ligation, reverse crosslinking/sub-library PCR) and
lines 3588–3597 (processing). The scraper (run on the article text) found the
Extended Data Fig. 1a structure only in OCR text of an image page; that OCR is not used
as a sequence source. No sequence or modification on the page changed.

- 🟢 Embryo nuclei: the same Homogenization Buffer, Dounce, 40 µm filter and NWB
  composition as GEO GSE243256; fixation by "adding 138 μl 37% formaldehyde (final
  concentration, 0.5%) for exactly 10 min on ice", 500 µl 2.5 M glycine, 5 min room
  temperature, two NWB washes, 5,000 nuclei per µl. Cell lines: Omni-ATAC lysis
  (0.1% NP-40, 0.1% Tween-20, 0.01% digitonin) 3 min on ice, then the embryo fixation.
- **Source discrepancy:** GEO GSE243256 prints 140 µl of 37% formaldehyde; the published
  Methods print 138 µl (0.5% final). The model uses the published value and the page
  shows both. Both give about 0.5% formaldehyde.
- 🟢 Transposome assembly: eight 20 µM 5′-phosphorylated Read 1 oligos and eight 20 µM
  5′-phosphorylated Read 2 oligos (each with an 8-nt barcode) annealed with an equal
  amount of 20 µM blocked ME-complement (5′P, 3′ddC), 85 °C 2 min, −1 °C/min to 20 °C;
  48 read 1/read 2 combinations chosen (Supplementary Table 1); 5 µl 10 µM transposon +
  10 µl Tn5 (Robustnique M0221) + 4 µl 10× TPS + 21 µl water, 37 °C 30 min; stored with
  glycerol at −20 °C. All oligos HPLC-purified (Sangon).
- 🟢 Ligation oligos: round-2 hybridization oligos share a universal linker; round-3
  oligos have 48 well-specific linkers. Linker and barcode strands (5 µl each, 20 µM) are
  pre-annealed in 10 mM Tris pH 8.0, 50 mM NaCl, 95 °C 2 min, −1 °C/min to 20 °C, giving
  a 5′ overhang complementary to the transposed end, a well barcode and a PCR handle.
- 🟢 Transposition: 10 µl nuclei (50,000) with MgCl2, 5× LM, NP-40, digitonin, protease
  inhibitor and 5 µl indexed Tn5, 37 °C 30 min at 300 rpm; stopped with 20 µl 40 mM EDTA
  + 1 mM spermidine, 37 °C 15 min. Species mixing used the first eight transposomes; the
  zebrafish atlas used all 48.
- 🟢 Rounds 2 and 3: pool, wash twice in NWB, resuspend in 1.06 ml ligation mix
  (1.5× T4 ligation buffer, 60 U/µl T4 DNA ligase NEB M0202L, 0.15% Triton X-100,
  0.36× NWB); 20 µl per well onto 10 µl barcodes; 37 °C 30 min at 50 rpm; stopped with
  EDTA/spermidine 37 °C 15 min. Species mixing used 10 U/µl ligase.
- 🟢 Sub-libraries: pool, wash, resuspend in Buffer EB, 1,000–5,000 nuclei per tube;
  proteinase K + 1% SDS 65 °C 16 h; 5% Triton X-100 quench 37 °C 1 h; PCR with i5 and i7
  index primers (Supplementary Table 1) in 2× KAPA master mix (KK2602): 72 °C 10 min,
  98 °C 45 s, 11 cycles (98 °C 20 s, 62 °C 30 s, 72 °C 30 s), 72 °C 5 min (13 cycles for
  species mixing); 1× SPRIselect cleanup; paired-end 150 bp on NovaSeq 6000 or DNBSEQ-T7.
- 🟢 Processing: genomic mates are raw Read 1 positions 66–121 and Read 2 positions
  51–101, as in GEO; rounds 1–3 barcodes are extracted from Read 1 and Read 2.
- 🔴 Not printed: conversion of the P5/P7 library for DNBSEQ-T7 runs.
- The split-pool model moved from the shared `lib/split_pool_atac.py` into
  `tools/spatac_model.py` (2026-10-17; SPATAC was its only caller). Its ligation steps now
  carry T4 DNA ligase NEB M0202L at {{= spatac_model.LIGASE_U_PER_UL }} U/µl in the
  ligation mix, {{= spatac_model.LIGATION_CELSIUS }} °C for
  {{= spatac_model.LIGATION_MINUTES }} min, and the PCR step carries the KK2602 master mix
  and {{= spatac_model.PCR_CYCLES }} cycles
  (species mixing: {{= spatac_model.SPECIES_MIXING_LIGASE_U_PER_UL }} U/µl ligase final,
  {{= spatac_model.SPECIES_MIXING_PCR_CYCLES }} cycles), each citing the Methods section.
  The nuclei-preparation steps cite the published Methods (with GEO as agreeing record).

## Second reading (independent, 2026-10-17)

Read directly in `s41556-024-01449-0.pdf.txt` lines 3364–3394 (nuclei isolation and
fixation), 3440–3559 (SPATAC-seq method through sub-library PCR), 3588–3597 (processing)
and the Extended Data Fig. 1 legend (line 4019), plus Extended Data Fig. 1a rendered from
PDF page 20 at 600 dpi. Supplementary Table 1 was re-checked for all 48 `Tn5_R_L` and 48
`Tn5_R_LA` rows.

- 🟢 Agrees with the first reading on transposome assembly, ligation-oligo annealing,
  transposition, the two hybridization–ligation rounds, reverse crosslinking, the PCR
  program and the 138 µl (Methods) vs 140 µl (GEO) formaldehyde discrepancy.
- 🟢 Extended Data Fig. 1a, read rendered: row 1 shows both round-1 adapters with a 5′
  phosphate (red star) on the leader ends; row 2 shows the round-2 oligo
  `5′-ACGACGCTCTTCCGATCT-[N8]-ATCCACGTGCTTGAG` joined to `AGGCCAGAGCATTCG` with the
  "left linker" `TAGGTGCACGAACTCTCCGGTCT` written 3′→5′ beneath it (the reverse of
  `Tn5_F_LA`); row 3 shows the "right linker" `CGCCGATG-[N8]-CTGTCTC` (`Tn5_R_LA-1`
  pattern) bridging the round-1 reverse leader and the round-3 oligo (S7 + ME + barcode,
  drawn 3′→5′); row 4 shows P5–i5–TruSeq Read 1 and P7–i7–S7 arms, a 65-nt Read 1 prefix
  and a 50-nt Read 2 prefix. These agree with the modelled construct.
- **Source discrepancy (figure-internal):** row 4 draws the round-3 barcode as
  `[NNNNNNN]` (seven N), while row 3, Supplementary Table 1 (all 48 `Tn5_R_L` oligos end
  in 8 barcode bases) and the figure's own 50-nt bracket require eight. The model uses 8;
  the page shows both values.
- 🟡 The ligation mix is 60 U/µl, and 20 µl of it is added to 10 µl of barcodes, so the
  final ligase concentration is about 40 U/µl (computed). The species-mixing value is
  printed as a final concentration (10 U/µl).
- Nuance (not harmonized): the Methods give the 66–121 / 51–101 genomic positions in the
  species-mixing processing paragraph; the zebrafish paragraph does not restate positions.
- The Methods name the PCR mix only as "2× KAPA Master Mix (KK2602; Roche)"; the model
  uses that wording and does not add a product name.

## Nuclei preparation (GEO sample protocols)

- 🟢 GSE243256 sample extraction protocol (zebrafish embryos): Dounce homogenization in
  Homogenization Buffer (10 mM Tris-HCl pH 7.8, 5 mM CaCl2, 3 mM Mg(Ac)2, 0.1% NP-40,
  0.1 mM EDTA, 320 mM sucrose, 0.015 mM PMSF, 0.15 mM β-mercaptoethanol), 40 µm filter,
  wash in nuclei wash buffer (NWB: 10 mM Tris-HCl pH 7.4, 10 mM NaCl, 3 mM MgCl2, 1% BSA,
  0.1% Tween-20); fixation by adding 140 µl of 37% formaldehyde to 10 ml NWB for exactly
  10 min on ice, quenched with 500 µl 2.5 M glycine for 5 min at room temperature; two NWB
  washes; nuclei diluted to 5,000 per µl in NWB. The published Methods print 138 µl
  (0.5% final) instead; see the discrepancy above.
- 🟢 GSE215901 sample extraction protocol (K562/Hepa): Omni-ATAC lysis (0.1% NP-40, 0.1%
  Tween-20, 0.01% digitonin, 10 mM NaCl, 3 mM MgCl2, 10 mM Tris-HCl pH 7.4; 3 min on ice),
  NWB wash, then "the fixation protocol of embryo nuclei above" (that protocol text is not
  in this series record), 5,000 nuclei per µl.
- Transposome assembly, transposition, ligation and PCR conditions are in the published
  Methods (section above), not in GEO.
- `tools/chem.py` `prepare_nuclei()` declares the embryo homogenization, formaldehyde
  fixation (published 138 µl value) and glycine quench, citing the published Methods with
  GEO GSE243256/GSE215901 as agreeing records, and passes them to the split-pool
  graph in `tools/spatac_model.py`.

## Bounded pass 2026-10-17 (round breadth-first-19)

Inspected: publisher HTML (subscription preview; abstract and Extended Data Fig. 1 legend,
lines 626 and 1529 of the text twin — the Methods are paywalled, manual URL
<https://doi.org/10.1038/s41556-024-01449-0>); Reporting Summary (MOESM2, image-only, rendered
p. 1 "Data collection"); Supplementary Table 1 (MOESM4 sheet 1, rows 4–52).

- 🟢 Supplementary Table 1 rows 4–14 and 52 match every representative oligo and
  modification on the page: `Tn5MErev` (5′P; 3′ddC), `Tn5_F_D1`/`Tn5_R_D1` (5′P),
  `Tn5_F_L1`, `Tn5_F_LA`, `Tn5_R_L1`, `Tn5_R_LA-1`, i5 primer-1 and i7 primer-1.
  No 5′ modification is printed for the round-2/round-3 oligos; their 3′ ends join the
  5′-phosphorylated round-1 adapters.
- 🟢 Reporting Summary: all libraries were sequenced paired-end 150 bp on Illumina NovaSeq
  6000 and HiSeq X Ten, or on MGI T7 and MGISEQ-2000. The 115/100 Read 1/Read 2 lengths on
  the page (set in the shared model) are the processing spans, not acquired cycle counts; the
  page now says so. Conversion of this P5/P7 library for MGI runs is not described in the
  inspected material.
- 🟢 Extended Data Fig. 1 legend: libraries use Illumina's standard recipe; the
  proof-of-concept used 8 barcoded transposomes and 8 PCR sublibraries.
- Open at that time: published Methods (then paywalled; read 2026-10-17, section above),
  GEO GSE215901 (reCAPTCHA; later read as SOFT), the MGI-conversion route (still open).

## Follow-up pass 2026-10-17 (round deferred-followup-02)

Routes tried for the then-paywalled Methods (later retrieved by hand through institutional access): Unpaywall and OpenAlex (closed; publisher and
PubMed locations only); Europe PMC (PMID 38977847, no PMCID, not in Europe PMC); Europe PMC
preprint search for SPATAC-seq and the zebrafish atlas (no preprint); the companion
mouse-organogenesis article doi:10.1038/s41556-024-01435-6 (also closed). Manual URL
remains <https://doi.org/10.1038/s41556-024-01449-0>. The GEO records, previously blocked by
reCAPTCHA on `acc.cgi`, were read as SOFT family files from the NCBI GEO FTP area (sections
above). The authors' repository file tree contains processing code, barcode lists and
figures only; no wet-lab protocol document.

## Bounded pass 2026-10-17 (round deferred-followup-03)

Hand-placed originals checked for identity: `s41556-024-01449-0.pdf` prints Nature Cell
Biology, doi:10.1038/s41556-024-01449-0, "Mapping the chromatin accessibility landscape of
zebrafish embryogenesis at single-cell resolution by SPATAC-seq", Sun, Liu, Xu … Lan
(received 15 Dec 2022, published 8 Jul 2024); MOESM1 prints the same DOI and title.
Re-read the Methods `s41556-024-01449-0.pdf.txt:3399-3559` (bulk ATAC, 10x scATAC, SPATAC
principle, transposome assembly, activity test, ligation-oligo preparation, transposition,
hybridization/ligation, reverse crosslinking and PCR) and processing lines 3588-3599.
They agree with the first and second readings above; no condition changed.

- 🟢 MOESM5 and MOESM6 (no text twins) are the Source Data workbooks: their sheet names
  are "Figure 1b" … "Figure 7k" and "Extended Data Fig.1f" … "Extended Data Fig.10h Right"
  (QC metrics, TSS enrichment, cell-type tables). They carry no oligo or reaction
  information. MOESM2 is the Reporting Summary (3 pages, image-only) and MOESM7 is
  Supplementary Figure 3, the uncropped gel of Extended Data Fig. 1c (1 page, image-only).
- 🟡 Model correction: the graph had discarded same-ended transposition products as
  products that "do not reach library PCR". Each well's transposome carries a read-1
  (`Tn5_F_D`) and a read-2 (`Tn5_R_D`) adapter, and the Methods describe no separation, so
  these fragments stay in the nuclei and enter the PCR tube. The discard reason now says
  that they are not physically separated and that their failure to amplify exponentially
  with one i5 and one i7 primer is computed, not observed.
- 🟢 The 10x scATAC comparison libraries (Chromium Next GEM ATAC v1.1, MGISEQ-2000,
  100 + 16 + 8 + 100) and the bulk ATAC libraries (Vazyme TD501, N5/N7 index primers,
  nine cycles, DNBSEQ-T7 PE150) are comparison assays, not SPATAC-seq, and are not drawn.
- Still open: the conversion of the P5/P7 library for DNBSEQ-T7 / MGISEQ-2000 runs.
