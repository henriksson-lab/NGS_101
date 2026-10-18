# DR-seq — research notes

## Sources read

- 🟢 Dey et al. 2015, DOI `10.1038/nbt.3129`, PMC author manuscript PMC4374170 (PDF text lines 68–92 main-text workflow; 330–435 Online Methods "DR-Seq" and sequencing; 439–470 bioinformatics) and the publisher Supplementary Information (MOESM36; Supplementary Fig. 3 legend lines 262–280; Supplementary Note lines 1358–1412).
- 🔴 The original CEL-Seq article and its supplement are not in the checked DR-seq source bundle. Teichlab's DR-seq reconstruction is useful for orientation, but it is not source evidence for Ad-1x bases or construct order.
- 🟢 Zong et al. MALBAC paper, DOI `10.1126/science.1229164`, for the looping/quasilinear mechanism.

## Shared amplification before the split

🟡 The working Ad-1x model is `CGATTGAGGCCGG` + T7 promoter `TAATACGACTCACTATAGGG` + RA5 `GTTCAGAGTTCTACAGTCCGACGATC` + an 8-nt cell barcode + dT24V. The DR-seq paper names Ad-1x and describes its cell-specific barcode and reverse transcription role, but it does not print this complete oligo. Its bases and length require verification in the cited original CEL-Seq material; the upstream reconstruction does not establish them.

🟡 The DR-seq paper prints `GTGAGTGATGGTTGAGGTAGTGTGGAG` as P2, which primes Ad-2-derived products. The main text (line 77) describes Ad-2 as a defined 27-nt sequence followed by 8 random nucleotides; the supplement (Supplementary Fig. 3 legend line 268; Supplementary Note line 1404, citing Zong et al.) says the 3′ end of Ad-2 is GGG or TTT. This published discrepancy is shown on the schematic. The working Ad-2 model keeps 8 variable positions, drawn as 5 random nt plus a named GGG-or-TTT terminal triplet; the complete Ad-2 primer mix is not printed in DR-seq.

🟢 Reagents (Methods lines 344–369): MessageAmp II first-strand reagents with ArrayScript, 42 °C, then 80 °C for 10 min; ThermoPol buffer and Ad-2; 94 °C 3 min; seven quasilinear cycles with fresh Bst large fragment plus Pyrophage 3173 exo- before each cycle; the 58 °C looping step is omitted in the first round. Seven quasilinear cycles randomly prime both genomic DNA and cDNA. Most short genome-derived products have Ad-2 at both ends; productive transcript-derived products have Ad-2 at one end and retain Ad-1x at the other.

🟢 The random priming coordinate—not a conventional UMI—is used as a length-based molecular identifier for RNA deduplication. The old page incorrectly presented the 8-nt Ad-1x cell barcode as an 8-nt UMI.

## Branch-specific finishing

🟢 The pool is split only after the shared amplification. The RNA aliquot receives P1 (`CGATTGAGGCCGGTAATAC`) for second-strand synthesis across the Ad-1x/T7 end in one PCR cycle; cells with non-overlapping barcodes are then pooled and column-purified before 13 hours of T7 IVT and aRNA clean-up without fragmentation. Genomic Ad-2/Ad-2 products lack that transcript-selective T7 architecture. Main text lines 89–90: 3′ Illumina adapters are ligated to the aRNA, followed by reverse transcription and PCR (details deferred to CEL-Seq, ref. 13).

🟢 The DNA aliquot is exponentially amplified with P2 (`GTGAGTGATGGTTGAGGTAGTGTGGAG`; Deep VentR exo-, 21 cycles), and 80 ng is re-amplified with 5′-biotin P3 (`GTGAGCTGGAGTTGAGGTAGTGTGGAG`; 4 cycles at 46 °C then 9 at 59 °C), sonicated to ~300 bp, and depleted of the biotin-bearing ends on MyOne C1 streptavidin beads (supernatant and two washes kept) before NEBNext Ultra library preparation with a different index per cell. P3 is the only primer listed for that PCR, so the model gives both product ends P3.

🟡 Printed P3 is not identical to P2: it differs at several internal positions near its 5′ end (computed on the schematic). The Methods do not explain this; both printed sequences are kept and the low first-cycle annealing temperature is recorded, without assuming a typo.

## Read layout

🟢 Methods lines 432–434: HiSeq 2500; DR-seq cDNA libraries 100 bp paired-end; gDNA and other cDNA libraries 50 bp or 100 bp paired-end (which length applies to which DR-seq gDNA library is not assigned). Lines 441–455: the left mate identifies the cell barcode; Ad-2 is trimmed from the right mate, which is mapped; the first right-mate coordinate after trimming is the length-based identifier.

🟡 The paper refers RNA library finishing to CEL-Seq and names NEBNext Ultra for DNA preparation. The cited CEL-Seq original, its supplement, and a version-matched NEBNext guide remain unreviewed. The full Ad-1x and Ad-2 oligos and final Illumina shells therefore remain provisional; only P1, P2 and P3 are exact source-printed primers in the checked DR-seq article.

## Open questions

- Complete Ad-1x and Ad-2 order sequences (CEL-Seq Hashimshony 2012 original and its supplement; MALBAC supplement not read).
- NEBNext Ultra adapter/index primers and which indices were used; gDNA acquired read length (50 or 100 bp).
- Whether the P3/P2 difference is deliberate.

## Native ownership migration, 2026-10-18

Read two originals: the published DR-seq article (PMC4374170 main workflow and Methods, PDF 69–93 and JATS 70–80) and its published supplement (Supplementary Figure 3 legend, publisher supplement 268–272; contextual published Ad-2 interpretation retained). Ran the cached-bundle oligo scraper and verified P1/P2/biotin-P3 against the original Methods. The cited CEL-Seq and MALBAC originals, historical kit manuals and source oligo orders were not followed in this capped pass.

Removed the shared multibranch Illumina endpoints and generic PCR finishers. Historical NEBNext Ultra DNA preparation is now an explicit bounded manufacturer operation; its actual end preparation, adapter/index and PCR recipe is not replaced with canonical TruSeq or a newer kit. The RNA branch follows the article's stated 3′-adapter ligation → library RT → PCR sequence, with exact historical CEL-Seq adapters/primers unassigned. Its PCR receives an actual single-strand library-cDNA Scene, distinct from the final library. DNA sample-index semantics retain the Methods' different index primers per single cell; the transcript barcode remains an Ad-1x cell identity, never a UMI.

Both halves of the quasilinear pool now retain genome- and transcript-derived populations. P2/P3 amplification does not imply physical modality purification; each PCR declares one physical primer species and uses its actual mixed incoming state, with first-cycle geometry unresolved. The initial RT product is an independent cDNA strand beside genomic DNA, rather than a primer/RNA docking view. Sonication produces a fragment population instead of an unchanged full P3 amplicon. Streptavidin depletion retains the source supernatant/wash DNA while excluding the biotin-labeled ends. DNA templates/nontranscribed carryover after IVT have unresolved recovery, rather than being presumed destroyed because only cDNA makes T7 aRNA. Previously labelled Ad-1x/Ad-2 working models and the published Ad-2 discrepancy remain explicitly inferred and are not promoted by this migration.

The RNA paired100 and DNA paired50-or100 HiSeq2500 alternatives are kept without assigning unreviewed index acquisition or exact native primer/read geometry. The source inventory is curator-written. No accepted missing-information claim or complete protocol review is implied.

### Independent component-inheritance correction

The later CEL-Seq library RT now has its own unassigned historical catalyst. ArrayScript/MessageAmp II is positively assigned to the initial cellular RT and is not reused as evidence for the later library RT. Both RT views explicitly leave the RNA template's subsequent physical disposition unassigned; a DNA product does not establish RNA destruction or removal. No new original was investigated for this implementation correction. The independent review otherwise found the mixed inputs, source physical fates, actual ss PCR precursor and vendor/geometry boundaries consistent within the stated source scope.
