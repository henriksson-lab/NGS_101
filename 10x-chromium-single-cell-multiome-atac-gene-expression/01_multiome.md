# 10x Chromium Single Cell Multiome ATAC + Gene Expression — research notes

## 1. What it is

The commercial Multiome assay jointly measures open chromatin and polyadenylated RNA from the same nucleus. Bulk nuclear tagmentation precedes GEM partitioning; the GEM barcode is then incorporated into both the ATAC and gene-expression products.

## 2. Sources read

- 🟢 10x Genomics, *Chromium Next GEM Single Cell Multiome ATAC + Gene Expression User Guide*, CG000338 Rev G (2024): [current support page](https://www.10xgenomics.com/support/instruments/chromium-x-series/chromium-next-gem-single-cell-multiome-atac-plus-gene-expression-reagent-kits-user-guide).
- 🟢 The guide's workflow, library diagrams and sequencing section establish bulk tagmentation, GEM barcoding, the two library branches, and their read roles.

## 3. Construct-changing steps

🟢 CG000338 Rev G Introduction (txt lines 1330–1400) and Steps 1–7 (lines 2100–4910):

1. Transposition in bulk nuclei: ATAC Buffer B + ATAC Enzyme B PN-2000265, 37 °C 60 min.
2. GEM generation on Chip J; Master Mix = Barcoding Reagent Mix, Template Switch Oligo
   PN-3000228, Reducing Agent B, Barcoding Enzyme Mix (reagent list lines 321–357, 2360–2390).
3. GEM incubation 37 °C 45 min, 25 °C 30 min, then Quenching Agent PN-2000269. One
   incubation produces barcoded cDNA (oligo-dT bead oligo + TSO) and 10x barcoded DNA from
   the transposed fragments (P5/barcode/spacer bead oligo). 🔴 The guide does not describe
   how the ATAC bead oligo is joined to the s5 end; the page draws a warned unresolved step.
4. Dynabeads MyOne SILANE then SPRIselect cleanup.
5. One pooled Pre-Amplification PCR (Pre-Amp Mix + Pre-Amp Primers PN-2000271; 72 °C 5 min,
   98 °C 3 min, 7 × [98 °C 20 s, 63 °C 30 s, 72 °C 1 min]) "to fill gaps" and amplify both
   modalities. 🔴 primer bases not printed.
6. Split: 40 µl to ATAC Sample Index PCR (SI-PCR Primer B + Single Index N Set A; 98/67/72 °C),
   35 µl to cDNA amplification (cDNA Primers PN-2000089; 63 °C anneal).
7. GEX: fragmentation + end repair + A-tailing, Adaptor Oligos PN-2000094 ligation, Dual Index
   TT Set A PCR (54 °C anneal). 🔴 adaptor strands and ends not printed.

Sequencing (lines 5145–5285): ATAC Read 1N 50, i7 8, i5 24 (8-nt spacer + 16-nt barcode; NextSeq
500/550 need a custom recipe with 8 dark cycles), Read 2N 49; GEX 28/10/10/90. 🟢

## 4. Evidence boundary

🟢 CG000338 Rev G Appendix pp. 80–81 prints the complete constant structure of both gel
bead oligos and both final libraries. The page now transcribes the GEX Read-1 arm,
16-base barcode, 12-base UMI, T30VN, dual 10-base sample indices, and the ATAC P5/barcode/
`CGCGTCTG`/s5/ME and s7/8-base-index/P7 shell exactly. Variable barcode, UMI and index
windows are structured placeholders, not inferred sequence. 🔴 The ATAC joining reagent
(if any splint) and the kit reaction-mixture oligos (transposon, TSO, Pre-Amp, cDNA,
adaptor and index primers) are not printed in Rev G.

### Guide-revision discrepancy

🟢 CG000338 Rev B Appendix p. 80 prints the gene-expression final library's 5′ P5 arm
as `AATGATACGGCGACCACCGAGATCT` before the `N10` i5 index. The ATAC gel-bead and
library lines on that same page include the following `ACAC`; Rev F and Rev G Appendix
pp. 80–81 also print the full `AATGATACGGCGACCACCGAGATCTACAC` in the
gene-expression library. 🟡 The model uses the full Rev G P5 arm, supported by the
other lines and canonical Illumina P5. 🔴 10x does not state whether Rev B's shorter
gene-expression line was a diagram error or represented a reagent variant.

## 5. Bounded source pass 2026-10-17 (round breadth-first-06)

Read CG000338 Rev G (reagent lists, Steps 1–7, Sequencing, Appendix pp. 80–81) and
CG000168 Rev A Appendix Step 4.1 for the SI-PCR Primer B part number. Model corrections:

- Graph order now follows Rev G: both GEM reactions, quench/cleanup and one pooled
  Pre-Amplification PCR precede the ATAC/cDNA split (the model had split before RT and
  had no pre-amplification, cDNA amplification, fragmentation or adaptor ligation).
- Template Switch Oligo, Pre-Amp Primers, cDNA Primers and Adaptor Oligos are declared by
  part number with unavailable bases/geometry; the ATAC bead join is a warned unresolved
  step rather than an extension by an undocumented "bridge oligo".
- The ATAC PCR geometry no longer claims that CG000338 prints the P5-side primer; SI-PCR
  Primer B is catalog-linked to its CG000168 Rev A sequence (same PN-2000128). 🟡
- Acquired cycles added: ATAC 50/8/24/49 (i5 reads spacer then barcode) and GEX 28/10/10/90. 🟢
- Same-adaptor tagmentation products and empty/multiplet GEMs now have unresolved fates
  instead of physical discard.

Still open: transposon and kit oligo bases; ATAC bead joining chemistry; identity of
Pre-Amp primer ends; CG000338 Rev B/F step differences beyond the Appendix comparison.

## Source-owned validation migration, 2026-10-18

Removed shared PCR-finishing recipe and its final-product primer-site assumptions. All four PCR geometries now refer to their actual preceding recovered/split/adapted inputs; the ATAC/GEX index PCRs remain explicitly unavailable because complete physical primers/precursors are unassigned, instead of claiming exact sites from the completed final library. The two local sequencing terminals retain source ATAC50/8/24/49 and GEX28/10/10/90 plans; the ATAC NextSeq500/5508dark+16reported cycle caveat is explicit. SI-PCRPrimerB is a source-qualified catalog-linked input with a direct predecessor-original transcription check. Same-end/GEM candidate classes have molecular unresolved fates, and the existing four-P5-base source-discrepancy notice remains. Known T30 is now a concrete segment separate from the two-base VN anchor, preserving sequence and exact lengths rather than treating all32positions as placeholders.

Source scope: 10x CG000338 RevG original: §5.1 SI-PCRPrimerB PN-2000128, §7.5 TT sample-index PCR, sequencing tables around lines5150–5178 and Appendix pp80–81 around lines5890–5950; catalog-linked CG000168 RevA Appendix SI-PCRPrimerB line3411. Existing RevB/RevG discrepancy source checks rerun without replacing their prior located reviews.

Deferred: Physical Pre-AmpPrimers PN-2000271 and cDNA/index primer mixtures, full transposome/adaptor states and ATAC barcode-joining mechanism remain unassigned; Complete intermediate PCR products and input-output molecular relations remain unresolved; RevB/RevF guide differences beyond the existing scoped P5 discrepancy were not re-reviewed; later revisions are not substitutes for cited RevG. This bounded check establishes neither complete chemistry validation nor confirmed source absence.
