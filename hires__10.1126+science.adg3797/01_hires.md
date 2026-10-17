# HiRES — research notes

## Sources read

- 🔴 Liu et al. 2023, DOI `10.1126/science.adg3797`. The Science full text and supplement are not in the source bundle (manifest `(manual)` row, https://doi.org/10.1126/science.adg3797); nothing on this page is taken from them.
- 🟢 Chen et al. 2023, DOI `10.21769/BioProtoc.4886`, [author detailed protocol](https://bio-protocol.org/pdf/bio-protocol4886.pdf), Materials p. 4 and procedure pp. 6–10 (extracted text lines 392–700: fixation C, RT D, digestion/ligation E, sort F, lysis/preamplification G, library H). Its abstract explicitly identifies the Science paper as the assay it implements.
- 🟢 Authors' GEO deposit GSE223917 (linked from PubMed 37289875): series record and sample GSM6998595 (`GEO_GSM6998595_self_brief.txt` lines 15–23, fetched by hand 2026-10-17). Not the Science methods, but author-written extract-protocol and data-processing text for the same study.

## Construct path and oligos

🟢 Step order and conditions (Bio-protocol procedure C–H, text lines 392–700): fix with 1.75% formaldehyde 10 min RT, BSA quench, Igepal lysis, 0.5% SDS 65 °C 10 min and Triton quench; in-situ RT with GAT5-RT and SuperScript IV (10 U/µl, 55 °C 15 min, betaine); MboI 37 °C overnight, 65 °C 20 min; T4 DNA ligase 16 °C 4 h in ~1 ml; FACS one nucleus per well of a dry 96-well plate; lyse in SC lysis buffer with Qiagen protease **and GAT5** (50 °C 3 h, 70 °C 1 h); preamplify with GAT5-7N + GAT5 and Deep Vent (exo-) (9 stepped 4→65 °C random-priming cycles, then 14 cycles at 58 °C annealing); tagment 2 µl with TTBL/TTE Mix V50 55 °C 10 min, stop with 0.2% SDS; 72 °C 3 min then 5 cycles with RNA-P7 + i5 UD index only ("slight enrichment of cDNA reads"); add i7 UD index and 5 more cycles; pool and select with 0.6× and 0.15× AMPure XP. The graph now declares the fixation, sort, protease lysis, the two-stage indexing PCR and the pooled size selection, which it previously omitted or merged.

🟢 The author protocol reverse-transcribes RNA **in situ before** MboI digestion and DNA proximity ligation. It then sorts individual nuclei into wells and jointly preamplifies contact DNA and transcript cDNA using GAT5-family primers. A Vazyme TruePrep V2 transposome fragments the mixed amplicons. RNA-P7 slightly enriches transcript reads before well-specific i5 and i7 dual indexes complete a **single joint library**. Contact and RNA reads are separated during analysis; they are not separate products of three cell-index ligation rounds. The previous page's three-round scheme was unsupported by this detailed source and has been replaced.

🟢 The four PAGE-purified primers printed in the Bio-protocol Materials section, with source line wraps rejoined, are:

- GAT5-RT: `5′-GTAGGTGTGAGTGATGGTTGAGGTAGTATTGCGCAATGNNNNNNNNTTTTTTTTTTTTTTTVN-3′`.
- GAT5-7N: `5′-GTAGGTGTGAGTGATGGTTGAGGTAGTNNNNNNN-3′`.
- GAT5: `5′-GTAGGTGTGAGTGATGGTTGAGGTAGT-3′`.
- RNA-P7: `5′-GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGGTTGAGGTAGTATTGCGCAATG-3′`.

🟢 GAT5-RT contains an eight-base random tract and anchored poly(T); GAT5-7N has a seven-base random 3′ tract. The protocol reports paired-end 150-base NovaSeq reads. The authors' [analysis repository](https://github.com/skelviper/HiRES) at commit `498941193bdb879c68ebb5d11f1157270ffcaf66`, `HiRES_preprocess_pipeline/rules/HiRES_split.rules` lines 1–37, classifies paired reads using the RNA-arm motif `XGGTTGAGGTAGTATTGCGCAATG` and then filters RNA pairs for the random-eight/poly(T) motif. This supports the two computational read classes from one joint library.

🟡 By computation, RNA-P7's 5′ portion is the reverse complement of the i7-side arm `CTGTCTCTTATACACATCTCCGAGCCCACGAGAC` in the TruePrep printed library structure, followed by the GAT5-RT arm `GGTTGAGGTAGTATTGCGCAATG`; it therefore primes on RNA-derived ends and supplies the arm an i7 primer can extend. Whether this matches the loaded Adapter 2 bases is not established.

🟢 The Vazyme TD501–503 manual (extracted PDF lines 145–178 and 295–311) says its tagmentation mix holds two adapter species at equal molar ratio and prints a standard finished TruePrep library with P5, i5, Read 1 arm, insert, Read 2 arm, i7 and P7. The [Illumina UD index reference](https://support-docs.illumina.com/SHARE/AdapterSequences/Content/Illumina-UDIndexes.htm) prints index adapter arms and plate-member sequences. 🟡 The schematic still marks the HiRES-specific shell inferred because the loaded transposome oligos, their exact terminal geometry, the RNA-P7 junction and chosen per-well i5/i7 members are not established by these sources. 🔴 The defining Science article and supplement remain unavailable in the source manifest.

## GEO sample-record text (round 2026-10-17-deferred-followup-02)

🟢 GSM6998595 extract protocol (lines 15–18): embryos dissociated in TrypLE Express, ~100 k cells in 200 µl PBS; cortex nuclei Dounce-homogenized in sucrose/KCl/MgCl2/HEPES/Igepal buffer with RNase and protease inhibitors, ~250 k nuclei in 500 µl Nuclei Storage Buffer. The library-stage paragraph (line 18) repeats the Bio-protocol: 15-fold dilution of preamplification product, 2 µl + 3 µl TTBL/TTE Mix V50 at 55 °C 10 min, 0.2% SDS stop; 72 °C 3 min and 5 cycles with i5 UD index + RNA-P7; add i7 UD index and 5 more cycles; pooled 0.6×/0.15× AMPure XP; paired-end 150-bp NovaSeq 6000. It prints RNA-P7 as `GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGGTTGAGGTAGTATTGCGCAATG`, identical to the Bio-protocol value (scraper hit, line 18).

🟢 Data processing (lines 19–21): RNA reads are those containing `GGTTGAGGTAGTATTGCGCAATG` "within the GAT5-RT primer from the 5-prime end of Read 2" (Cutadapt 3.5); reads with no poly-A after the GAT5-RT sequence are removed; Read 1 is aligned with STAR; duplicates are removed "based on UMIs with UMI-tools". 🟡 This supports, but does not by itself define base-for-base, the 8-nt GAT5-RT random tract as the UMI.

Routes tried for the Science article and supplement, all unsuccessful: OpenAlex/Unpaywall (closed, no repository copy); Science supplement URL `https://www.science.org/doi/suppl/10.1126/science.adg3797/suppl_file/science.adg3797_sm.pdf` HTTP 403 (not circumvented); Europe PMC search (no PMC copy, no preprint record). Manual URL remains https://doi.org/10.1126/science.adg3797. Still deferred: Science methods/supplement, TruePrep loaded adapter bases, selected UD index members and index cycles.
