# GRO-seq — research notes

Core, Waterfall and Lis, *Science* (2008),
[10.1126/science.1162228](https://doi.org/10.1126/science.1162228).

## Originals read

- 🟢 Original article, PMC2833333, particularly main text lines 28–32: BrUTP run-on,
  triple antibody selection and two biological replicates sequenced as 33-base reads.
- 🟢 The **original 2008 Supporting Online Material has been recovered** from a public
  mirror of the original PDF. Its cover names the defining title, authors, publication
  date, DOI and Science DC1 URL. Source file: `594009_9144e1ea58ad4573871ec024827290f3.pdf`.
- 🟢 Original SOM overview pp.9–10, Methods pp.17–19 (text 629–714), read processing
  (715–746), Figure S1 caption and original full figure on PDF p.27, and S15/S16 captions.
- 🟢 Later detailed protocol PMC5880289: actual oligo list and poly(A)-priming/circularization
  method were read as a **different implementation**. Its NTI223/200/201/202 oligos are
  not assigned to the original 2008 library.

2026-10-17 bounded re-reading (round 2026-10-17-breadth-first-11): SOM overview lines
381–410 and Methods lines 629–714 re-read against the local graph; article JATS lines
28–32 re-checked for the read layout; the scraper over the whole bundle returned four
hits, all in the later PMC5880289 implementation (NTI223/200/201/202) and none in the
2008 originals, confirming that no 2008 oligo bases are printed in the available text.

Mandatory oligo scraping preceded source reading. The original supplement's former
“unavailable” record describes a gated PMC route, not the availability of every original
copy. That history is retained in the source ledger; it no longer blocks reading this SOM.

## Original 2008 reaction order

1. 🟢 Isolate IMR90 nuclei by cold swelling, Igepal lysis and washing.
2. 🟢 Already-engaged polymerases extend nascent RNA with ATP, GTP, BrUTP and limiting
   CTP in Sarkosyl-containing buffer. The original includes α-32P-CTP as a purification
   tracer. BrU and the tracer are chemical labels, not molecular identifiers.
3. 🟢 DNase I digestion is followed by SDS/proteinase K treatment, acid phenol/chloroform
   extraction and ethanol precipitation of RNA.
4. 🟢 Base-hydrolyse RNA, neutralize it, pass it twice through P-30 columns, then perform
   a **second DNase I digestion** and EDTA stop.
5. 🟢 First anti-deoxyBrU capture: binding/washes, collection of four eluates, extraction
   and precipitation recover the BrU-containing RNA fraction.
6. 🟢 End repair is sequential: TAP removes remaining caps; low-pH PNK removes the
   hydrolysis-derived 3′ phosphate; raised-pH PNK with ATP adds 5′ phosphate. Extract
   and precipitate the repaired RNA.
7. 🟢 Ligate the **5′ adapter first**, using T4 RNA ligase-1 and PEG, then perform the
   second antibody capture/elution/extraction.
8. 🟢 Ligate the **3′ adapter second**, then perform the third antibody capture/recovery.
9. 🟢 SuperScript III reverse transcription uses the historical kit RT primer in two
   separate 10-µL reactions. RNase cocktail plus RNase H then digest RNA before PCR.
10. 🟢 Phusion amplification uses 15 cycles and the Illumina-specified PCR primers.
    Native 8% TBE PAGE selects cDNA greater than 90 nucleotides; elution, extraction
    and precipitation recover the sequencing library (SOM Methods lines 703–709).

🟢 SOM overview lines 403–410: the finished library averages **about 100 bases**
(“~190 base – 90 base adapters”), each of the two libraries had a pre-amplification
complexity of 1 × 10⁹ molecules estimated by re-amplifying a known aliquot, and 50
molecules were conventionally cloned and sequenced as quality control before the
Illumina 1G Genome Analyzer run. The ~90-base adapter contribution **bounds** the
combined length of the unassigned historical 5′/3′ adapter and primer bases; it does
not identify them, so the model still leaves those bases unassigned.

🟢 SOM Methods lines 654–657: the run-on reaction buffer supplies 500 µM ATP, GTP and
BrUTP with 2 µM CTP plus 0.33 µM α-32P-CTP, mixed with an equal volume of nuclei and run
5 min at 30 °C. The overview's “CTP is adjusted to 1 µM” (line 387) is the same value
after the equal-volume mix, not a conflicting figure.

The local graph owns this route and all three captures. Removed binding/wash, extraction,
column and size-selection fractions are represented separately. Digestion products are
not automatically called discarded material at the enzyme step. The repeat capture helper
implements only this protocol's own repeated source recipe.

## Historical adapter boundary

The SOM names **Illumina Small RNA Isolation Kit** 5′/3′ adapters, RT primer and PCR primers.
The current model leaves exact ordered bases, modifications, kit revision and primer sites
unassigned pending historical source reconciliation. This is an open model question, not
an accepted claim that these details were never published.

The later PMC5880289 method uses poly(A) tailing, oligo(dT)-containing NTI223, circularization
and abasic-spacer cleavage. Those operations and oligos are not a reconstruction of the
2008 ligation route. The earlier repository's generic “RNA adapter → RT → PCR” graph and
44-base insert/eight-base “UMI” label were modelling shortcuts; their removal is a repository
correction, not a publication error.

## Sequencing and computation

🟢 The original Illumina 1G experiment produced **33-base reads** from the 5′ end of
fragmented BrU-containing RNA. The read 5′ coordinate was used for polymerase-position
analysis. It is not an exact identification of a known current sequencing-primer family.

🟢 Iterative computational trimming from the read 3′ end down to a minimum of 29 bases
increased unique alignments. rRNA alignments were removed computationally. These are
post-sequencing analysis operations and do not physically trim or discard library molecules.

## Remaining original-source scope

The source checklist remains in progress. Actual obligations include the article-linked
GSE13518 metadata, nuclei-isolation predecessor Strobl and Eick (1992; SOM reference 39),
PNK 3′-phosphatase provenance Cameron and Uhlenbeck (1977; reference 21), and identification
of the historical kit guide. Full result-figure/table and remaining linked-reference scope
is not claimed from the bounded original Methods/Figure S1 reading. No missing-information
claim is accepted at this checkpoint.
