# MATQ-seq — research notes

> 🟢 source transcription · 🟡 derived · 🔴 unavailable

Defining source: Sheng K et al. “Effective detection of variation in single-cell
transcriptomes using MATQ-seq.” *Nature Methods* 14, 267–270 (2017).
[doi:10.1038/nmeth.4145](https://doi.org/10.1038/nmeth.4145)

Sources checked in the local bundle and linked originals:

- 🟢 the complete supplementary MATQ-seq protocol (`...MOESM227...pdf.txt`);
- 🟢 the [author-uploaded facsimile of the defining article](https://www.researchgate.net/publication/312482780_Effective_detection_of_variation_in_single-cell_transcriptomes_using_MATQ-seq), Online Methods, especially “Reverse transcription and second-strand synthesis,” “Library preparation and sequencing,” and “Overall cost per cell”; it specifies the adapter family, end preparation and run format absent from the local supplementary protocol;
- 🟢 Zong and Sheng's [2017 protocol preprint](https://doi.org/10.1038/protex.2016.088), pp. 1–3, which independently prints the same six oligos and reaction sequence; its original PDF is listed in `catalogue/extra_sources.tsv`;
- 🟢 Illumina's [Adapter Sequences v00 (October 2015)](https://dnatech.ucdavis.edu/sites/g/files/dgvnsk15531/files/inline-files/illumina-adapter-sequences_1000000002694-00.pdf), pp. 15–17, which lists separate TruSeq HT dual-index and LT single-index adapter families; the article does not select one;
- 🟢 supplementary figures, especially the paired-end insert/read description;
- 🟢 the individual publisher supplementary data workbooks MOESM217–219 and MOESM221–225, plus MOESM226 MATLAB analysis archive, inspected for chemistry relevance. These contain expression/count/analysis results rather than oligo or reaction instructions. The source review records the sheet and archive inventory.

## Reaction order

🟢 Single-cell lysate contains GAT27dT, GAT27-5N3G and GAT27-5N3T. SuperScript III
runs ten stepped-temperature cycles (8–50 °C), allowing repeated random and poly(A)
annealing. T4 DNA polymerase digests unused primers, optional Exonuclease I follows,
and RNase H and RNase If remove RNA (defining Online Methods, lines 1305–1327).

🟢 TdT adds dC to first-strand cDNA. `GAT21-6N3G` uses terminal GGG against that tail,
followed by ten Deep Vent exo− annealing/extension cycles. A single GAT27 primer then
amplifies for 24 cycles. `3NGAT24` converts the amplified material to more completely
double-stranded product before shearing.

🟢 All six protocol oligos are transcribed verbatim into `tools/matq_seq.py` and rendered
as raw selectable ordering text.

## Final library boundary

🟢 The defining article's Online Methods (lines 1369–1404) specify Covaris shearing to
approximately 200 bp, NEBNext end repair, Klenow exo− dA tailing, ligation of Illumina
TruSeq adapters with T4 DNA Quick Ligase, five cycles of KAPA HiFi PCR with a TruSeq PCR
primer, and NextSeq 500 sequencing. The exact selected TruSeq adapter/index oligo and
ligation-ready end geometry are 🔴 unavailable in the checked materials. The page
now uses source-limited TruSeq endpoint roles with no canonical terminal bases or invented index widths; native library-PCR geometry is unavailable on its actual adapter-ligated precursor.
The author protocol preprint links a `supplement0.docx` attachment that remains
inaccessible; the source review keeps this missing-detail claim blocked.

**Published read-length discrepancy.** The defining article's Online Methods, “Overall
cost per cell” (lines 1475–1485), explicitly say a NextSeq 500 **2×75** midoutput kit.
Supplementary Figure 21's caption (MOESM220 text lines 312–319) instead calculates
library fragment size using **85×2** bases of paired-end reads. The page models the
named 2×75 kit. The reason the supplementary caption says 85×2 and whether it describes
a different run are unresolved.

## Bounded curation round 3 (2026-10-16)

Inspected: Published Supplementary detailed protocol MOESM227 PDF text1–120; MOESM217 XLSX heading identifies coverage statistics.

Distinguished24-cycle WTA and subsequent20-cycle3NGAT24 conversion in the summarized graph label. Source confirms T4DNApolymerase primer digestion, RNaseH/I and TdT dC tailing before second-strand synthesis.

Deferred: Published main article and primer table, complete adapter/index acquisition and source mixture-unit wording; graph separation of WTA/conversion and recovery. This capped pass does not establish a complete source review or publication-wide absence.

## Native graph migration (2026-10-18, batch02)

Inspected: Published supplementary protocol MOESM227 extracted text1–142 (entire protocol); published supplementary figures MOESM220 Fig21 caption309–319. Mandatory whole-bundle scraper reviewed; six primer sequences verified directly in the original supplementary protocol.

Native graph separates three RT primer/target populations, retained RNA, source RNase digestion, variable TdT dC tail, ten-cycle second-strand extension, 24-cycle GAT27 WTA and two 20-cycle 3NGAT24 conversion blocks with their separate recoveries. Fragmentation, end repair, dA tailing, source-limited TruSeq ligation and five-cycle KAPA PCR use actual precursor states. Prior independently reviewed Online Methods end preparation and five-cycle/2×75 run scope retained, not re-read in this pass. The earlier bounded conversion label understated the two source blocks. Fig21 85×2 versus named 2×75 run discrepancy remains displayed.

Deferred: Full adapter/index family and ligation geometry; complete WTA/conversion molecular geometry; cached defining main article not re-read this pass; acquired index-cycle assignment. This capped pass does not establish publication-wide absence or complete molecular proof.
