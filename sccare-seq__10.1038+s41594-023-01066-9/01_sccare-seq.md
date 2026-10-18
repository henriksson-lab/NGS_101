# scCARE-seq — source notes

Defining paper: Qu J, Sun J, Zhao C, … Fan L, Ding J. Simultaneous profiling of chromatin architecture and transcription in single cells. *Nature Structural & Molecular Biology* 30, 1393–1402 (2023), [doi:10.1038/s41594-023-01066-9](https://doi.org/10.1038/s41594-023-01066-9). (Earlier versions of this note and page wrongly named "Deng et al."; corrected 2026-10-17.)

- 🟢 Defining article read (deferred follow-up 03, 2026-10-17). `s41594-023-01066-9.pdf`
  (publisher PDF, hand-downloaded; printed DOI, title and authors checked against the
  catalogue DOI). Text twin `s41594-023-01066-9.pdf.txt`: Results lines 232–250; Methods
  "Tn5 preparation" lines 2335–2370 and "scCARE-seq experimental procedure" lines
  2372–2465; "CARE-seq" lines 2467–2485; data preprocessing lines 2533–2545; Extended
  Data Fig. 1 caption line 2757 (schematic itself not inspected as image).
- 🟢 Methods step order and conditions (lines 2372–2465): 1% formaldehyde 10 min, 125 mM
  glycine; IGEPAL lysis to nuclei; RT of 1.6 × 10^4 nuclei with 8 µl 25 µM
  oligo-dT primer "5′-ACGACGCTCTTCCGATCTNNNNNNNN[8 bp RNA barcode]T30VN-3′", SuperScript IV,
  4→55 °C 2-min steps then 55 °C 10 min; 0.6% SDS 62 °C 10 min, Triton X-100 quench,
  MboI (25 U/µl) in NEBuffer 2, 37 °C 3 h; ligation with T4 DNA ligase (M0202S),
  16 °C 4 h; **second-strand synthesis** (NEB E6111) 16 °C 180 min; 40 µm strainer, DAPI,
  FACS one nucleus per well; lysis in 3 µl scReverse Crosslink Buffer (proteinase K
  1 mg/ml, **0.5 µM carrier ssDNA**) 55 °C 2 h / 65 °C 12 h / 70 °C 30 min; tagmentation
  with 0.025 µM in-house Tn5 in TAPS–DMF 55 °C 10 min, 72 °C 5 min, then proteinase K;
  pre-amplification with a mix of common, DNA and RNA primers (40 µM each) in Q5 2×,
  72 °C 5 min then 12–14 cycles; ExoI 37 °C 60 min / 80 °C 20 min; split 12 µl + 12 µl;
  Hi-C library: Nextera i5 + Nextera i7, 4 cycles (60 °C anneal); mRNA library: TruSeq i5 +
  Nextera i7, 10 cycles (65 °C anneal); each pooled, Zymo DCC-5 then 0.8× AMPure XP;
  "sequenced on the Illumina platform" (no read lengths).
- 🟢 Bulk CARE-seq only (lines 2467–2485): oligo-biotin-dT
  "5′-ACGACGC[BiodT]CTTCCGATCT…", Klenow fill-in with biotin-14-dATP, tagmentation without
  sorting, biotin enrichment. The single-cell procedure lists none of these, so the model
  has no fill-in, no biotin and no streptavidin step.
- 🟡 Contact junction: MboI ends ligated without fill-in reconstitute a single GATC
  (derived from enzyme specificity; the article prints no junction).
- 🟡 Pre-amplification pairing: N5-PreA/N7-PreA are the Tn5 S5/S7+ME 5′ portions and
  A5-PreA contains the RT handle, so contact fragments are inferred to be copied by
  N5-PreA/N7-PreA and cDNA by A5-PreA/N7-PreA; the products are not drawn base by base.
- 🟢 Extended Data Fig. 1 caption confirms "DNA/RNA libraries construction was worked by
  specific DNA or RNA primers"; the branch-specific index primers now come from the
  Methods text (Nextera i5 for Hi-C; TruSeq i5 for mRNA).
- Model corrections (2026-10-17, follow-up 03): authors corrected to Qu, Sun et al.;
  second-strand synthesis added before sorting; lysis/reverse-crosslinking with the
  carrier ssDNA made its own step; ExoI primer digestion added after pre-amplification
  (degraded primers left in the well, not discarded); contact junction drawn as GATC
  instead of an unspecified bridge; Oligo-dT annotation no longer suggests biotin in the
  single-cell protocol; **A5-PreA transcription corrected** from `CCCTAACGACGC…` (a missing
  C, repository error found by the source-oligo guard) to the printed
  `CCCTACACGACGCTCTTCCGATCT`; pre-amplification and index primers carry exact states.
- 🔴 (unread, not confirmed absent) Open: the role of the LIANTI carrier ssDNA; whether the
  Tn5 is loaded with a barcoded adapter (Tn5 preparation text says "annealed barcoded
  adapter", Supplementary Table 1 prints unbarcoded Tn5-N5/N7; model uses the printed
  adapters); read lengths and acquired index cycles; Protocol Exchange / predecessor
  methods (refs 34, 46, 55, 56, LIANTI).
- 🟢 Supplementary Table 1 prints Tn5-N5, Tn5-N7, phosphorylated Tn5-Rev,
  carrier ssDNA, the common/DNA/RNA pre-amplification primers, all eight
  UMI–barcode–oligo(dT) RT primers, and the complete N501–N516, N701–N724 and
  T501–T516 indexed PCR-primer sets.
- 🟢 Supplementary Table 1 row 21: "T: if it is Biotin dT, named Oligo-biotin-dT".
  The number and position of biotin-dT residues are not stated there; GEO GSE211395 names oligo-biotin-dT with biotin enrichment for bulk CARE-seq only
  (see below); the single-cell model no longer draws a streptavidin separation.
- 🟢 Peer Review File (MOESM2): reviewer #1 (lines 139–141, 151–152) describes
  "cell fixation [1% formaldehyde] followed by sequential reverse transcription and
  Hi-C"; reviewer #2 (lines 182–183) reads Fig. 1a as fixation of the whole nucleus,
  then RT; the authors' responses (lines ~395–400, ~480–483) say they captured the
  3′ end of mRNA to mitigate fixation effects. Reviewer text is not the authors'
  Methods; GEO (below) now gives the authors' own order.
- 🟢 GEO GSE211395 sample protocols (deferred follow-up 2026-10-17). The series links PubMed 37580628. The authors' deposited SOFT text is cached as `geo_GSE211395_acc-targ-gsm_brief.soft.txt`. Single-cell samples, for example GSM6468394 (`2i_H_1`), have `extract_protocol` at lines 447–448: "Cells were fixed with 1% formaldehyde for 10 minutes. Cells were flow sorted into individual cells before cell lysis, and chromatin was tagmentation with Tn5." and "mRNA reverse transcription, genome digestion (MboI) and ligation were performed in intact nuclei, and the individual cells were sorted by flow cytometry into separate wells and lysed. After adding Tn5 transposase to each well for tagmentation, each cell library was pre-amplified, and divided into two halves for DNA and RNA library amplification using respective primer pairs." The instrument is NovaSeq 6000 (line 466). Bulk CARE-seq samples (`CARE_rep1`, lines 272–273) say "oligo-biotin-dT was used for reverse transcription, and Tn5 tagmentation was performed after cell lysis without single cell sorting. The libraries were constructed after biotin enrichment."
- Model correction from that text. The order is now RT → MboI digestion → proximity ligation (all in intact nuclei) → flow sorting into wells and lysis → Tn5 tagmentation per well → pre-amplification → split into two halves → DNA (N5xx/N7xx) and RNA (T5xx/N7xx) indexed PCRs. Tagmentation after ligation is now source-stated, not inferred.
- Removed: the earlier streptavidin separation of a biotin-cDNA branch. GEO describes biotin enrichment only for bulk CARE-seq; the single-cell branches are halves of one pre-amplified library.
- 🟢 The printed N7-PreA/N5-PreA/A5-PreA primers (Supplementary Table 1 rows 8–10) are the Methods' "common primer, DNA primer, RNA primer" pre-amplification mix.
- 🟡 The eight RT barcodes are reclassified from cell barcode to a per-RT-reaction sample index, because RT precedes single-cell sorting (GEO). The model infers that per-cell identity comes from each well's i5/i7 pair.
- 🟢 Methods "Hi-C library preparation" / "mRNA library preparation" assign Nextera i5 (N5xx) + N7xx to the DNA-contact branch and
  TruSeq i5 (T5xx) + N7xx to the RNA branch. Thus the RNA endpoint is a hybrid library: TruSeq Read 1
  and i5 geometry at one end, Nextera Read 2 and i7 geometry at the other.
- 🟡 The page shows Oligo-dT-1, N501, N701 and T501 as exact representatives;
  the remaining printed members change only the documented whitelist sequence.

## Review trail

2026-10-17 deferred follow-up (round deferred-followup-02): no open copy of the defining article found (OpenAlex W4385812298 and Unpaywall report it closed, no PMC/Europe PMC full text, no preprint found in a Europe PMC search). Read the GEO GSE211395 sample protocols and modelled them as above.

2026-10-17 deferred follow-up 03: defining NSMB article now on disk (hand-downloaded publisher PDF) and read as listed above; corrections recorded above. Historical review trail entries before this one describe the state when the article was unavailable.


## Critical re-review, 2026-10-17

Candidate tagmentation products lacking the selected library ends remain unresolved; no represented separation establishes physical discard or amplification behavior.

Inspected scope: Local tagmentation-to-PCR declarations and terminal-fate audit; no new reaction-source assertion. This bounded correction does not complete the linked-source or two-reader absence audit.


## Native graph migration, 2026-10-18 batch05

Re-read two published originals: Qu, Sun etal.2023 main PDF, Tn5 preparation
and “scCARE-seq experimental procedure” (text2335–2465, fixation through
quantification); Supplementary Table1 MOESM3 workbook, oligo/index rows1–84.
The full bundle primer scraper was run, and hits checked against the table
including the phosphorylated Tn5-Rev, plain single-cell RT oligo, carrier and
three preamplification primers. Bulk biotin-dT is a separate route.

🟢 Local graph operations now explicitly retain RT→MboI→T4 proximity ligation→
NEB second strand→FACS→proteinaseK reverse-crosslinking with independent carrier
ssDNA→Tn5→three-primer Q5 preamplification→ExoI→two mixed aliquots→source
index PCR→pooled DCC-5 column→0.8×AMPure. The DNA and RNA aliquots both
contain the DNA/cDNA/carrier mixture. Different primer sets establish intended
target amplification, not physical modality purification.

🟡 Scene populations distinguish ss RT cDNA, complementary second-strand cDNA,
chromatin/contact DNA, native RNA and independent carrier. Native RNA and carrier
recovery/amplification remain unassigned. Actual incoming tagged/preamp states
are used for PCRGeometry.unavailable; final source-primer-derived libraries are
separate references, without claiming native first-cycle/gap/end proof. ExoI
hydrolysis remains chemically distinct from later physical column/bead removal.

The published Tn5 “barcoded adapter” wording versus unbarcoded Table1 oligos
remains visible. Run cycles and acquired index roles remain unassigned; exact
N8 physical index primers do not establish acquisition. This two-original cap
does not finish historical linked-source or absence review.


Implementation correction after independent review: tagged cDNA now explicitly retains a Tn5-derived terminal role through the tagged/preamp mixed populations. Transfer extent and gap geometry remain unassigned; no new bases or physical purification are inferred.
