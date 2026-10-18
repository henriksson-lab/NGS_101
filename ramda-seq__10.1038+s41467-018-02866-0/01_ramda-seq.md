# RamDA-seq — source notes

Defining paper: Hayashi et al., *Nature Communications* 2018, [doi:10.1038/s41467-018-02866-0](https://doi.org/10.1038/s41467-018-02866-0).

## Sources read

- 🟢 PMC5809388 full-text XML, Methods “RT-RamDA”, “Not-so-random primers”, “RamDA-seq”, and “Sequencing library information”.
- 🟢 Supplementary Note 1 and Supplementary Figure 1 in the paper's main supplementary PDF.
- 🟢 Supplementary Data 3 (`s41467-018-02866-0_supp_41467_2018_2866_MOESM6_ESM.xlsx`),
  the complete first- and second-strand mouse NSR ordering lists.

## Evidence trail

- 🟢 The first-strand pool contains 408 individually synthesized mouse NSR hexamers chosen to exclude exact rRNA matches; oligo(dT)18 is also present.
- 🟢 A representative matched pair is `1st_MmNSR_1` = `AAAGCA` and
  `2nd_MmNSR_1` = `TGCTTT`; the latter is the reverse complement in the published
  Supplementary Data 3 table. The page shows this pair as an ordering example,
  not as a replacement for the 408-member pools.
- 🟢 DNase I nicks cDNA in the RNA/cDNA hybrid. RNase-H-minus PrimeScript extends from nicked 3′ ends and displaces downstream cDNA; T4 gene 32 protein promotes displacement and protects displaced products.
- 🟢 Complementary second-strand NSRs and Klenow Fragment (3′→5′ exo−) make double-stranded cDNA.
- 🟢 Scaled Nextera XT tagmentation and 13 PCR cycles produce the original library; the advanced method uses 14 cycles.
- 🟢 Defining datasets used 76-cycle single-read NextSeq 500 sequencing. The protocol has neither cellular barcodes nor UMIs; the authors explicitly identify strand-displacement amplification as preventing their addition.

The 408+408 individual NSR sequences are supplied in Supplementary Data 3; the
full pool is not reproduced as an ordering list here.

## Bounded curation round 4 (2026-10-16)

Inspected: Published PMC5809388 JATS180–204, standard/advanced/C1 RamDA preparation and sequencing.

Added initial DNase I genomic-DNA digestion before RT and explicit unresolved fragments. Specified 13-cycle standard versus 14-cycle advanced library PCR; candidate Tn5 products no longer presumed discarded. Confirmed single-read76 acquisition.

Deferred: Original NSR tables, pre-/post-tagmentation cleanup, concurrent nicking/RT representation; C1 variant is separate and its JATS 130°C wording needs PDF verification before interpretation. This pass does not establish a complete source review or publication-wide absence.

## Native graph migration (2026-10-18, batch04)

Inspected three published originals: PMC5809388 JATS Methods175–204 (standard,
advanced, C1 context and sequencing), original MOESM1 Supplementary Note1 pp.1–2
(coupled RT/nicking/displacement), and original MOESM6 Supplementary Data3
`Mouse NSR` rows1/410 (representative complementary ordering pair). The bundle
scraper returned251 hits in3hit-files; hits are reading pointers, not facts.
The full408-member pool was not independently transcribed in this pass.

The graph selects the standard plate13cycle method. DNaseI digestion and ERCC
spike-in precede a coupled RT-RamDA reaction with retainedDNaseI, PrimeScript and
gp32, not separate reagent additions for completedRT then nicking. SampleRNA,
ERCCRNA, nickedDNA and displacedssDNA are independent populations. Klenow makes
a role-level duplex beside carriedRNA/fragments. Source15µlAMPure recovery and
direct6µlTagmentDNA Buffer elution precede scaledXT tagmentation. Actual native
transferred strand/gap classes enter13cyclePCR, then1.2×AMPure/24µlTE recovery.
Nonamplifiable products have unresolved physical fates. IndividualRNA/hydrolysis
partition through beads is unassigned; the retained cDNA drawing is not a claim
of complete removal of those species.

The old canonical XT shell/index widths and exact primer set were model assumptions.
Native ends/index features now have unknown encoding/width; positive single-read76
is retained without guessed index acquisition. Source-linked exact historicalXT
manual/loadedstrands, completeNSRpools, advanced andC1 reconstruction remain deferred.
No originalPDF resolution of C1's printed130°C wording was attempted. Two prior
constant-equals-itself tests were replaced by source transcription checks for the
representative complementaryNSRs; missing external tables skip those tests.
