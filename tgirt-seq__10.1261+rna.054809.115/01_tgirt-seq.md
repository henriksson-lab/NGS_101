# TGIRT-seq — original Qin 2016 ligation method

Defining source: Qin *et al.*, *RNA* 2016,
[doi:10.1261/rna.054809.115](https://doi.org/10.1261/rna.054809.115).
The diagram represents this original method, using its Supplemental Figure S1.

## Original oligos and molecular order

- 🟢 Figure S1 prints the 34-nt blocked R2 RNA and complementary 35-nt R2R DNA
  with an equimolar single-base DNA overhang. The overhang pairs to the target
  RNA’s terminal base; it becomes part of the target-derived cDNA, not a UMI.
- 🟢 The source RNA is
  `rArGrArUrCrGrGrArArGrArGrCrArCrArCrGrUrCrUrGrArArCrUrCrCrArGrUrCrArC/3’SpC3/`.
  The figure typesets the spacer as `/3’SpC3/`; Methods names the same C3 spacer
  `3SpC3`, used as the model’s chemical designation.
- 🟢 R2R is `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTN`.
- 🟢 Ordered R1R is
  `/5Phos/GATCGTCGGACTGTAGAACTCTGAACGTGTAG/3SpC3/`.
- 🟢 The original unindexed P5 PCR primer is
  `AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGACGATC`.
  The original P7 primer is
  `CAAGCAGAAGACGGCATACGAGAT[BARCODE]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`.
  The diagram leaves the selected barcode and its width unassigned.
- 🟢 Methods, “Construction of RNA-seq libraries” (JATS line 191), describe
  NaOH RNA hydrolysis, neutralization, cDNA cleanup, ligation of **5′-adenylated,
  3′-blocked** R1R to the cDNA 3′ end, and Phusion PCR. The ligated strand is
  therefore **5′ R2R — cDNA — R1R 3′**. The source P5 and P7 annealing arms face
  inward on that pre-PCR template, which the model checks directly.
- 🟡 Figure S1’s ordered 5′-phosphorylated R1R and Methods’ activated 5′-adenylated
  input are separate states. The page exposes an unresolved activation step and
  does not assign a particular adenylation reagent or enzyme to 2016. This is a
  model boundary, not a demonstrated publication error.

## Enzymes and sequencing scope

🟢 The original Methods distinguish TeI4c-ΔEn, TeI4c-MRF and GsI-IIC/TGIRT-III
across datasets (JATS line 175 and Supplemental Table S4). The graph keeps that
source variation explicit rather than assigning TGIRT-III to every experiment.

🟢 Plasma libraries were sequenced with 100-nt paired HiSeq 2500 reads or
75-/150-nt paired NextSeq 500 reads (Methods line 191). The page states these
alternatives. It does not assign one read length to all datasets or infer the
loaded sequencing-primer oligos from PCR-binding sites. Index-2 acquisition is
unassigned; absence of an i5 window in the original P5 primer does not by itself
prove which technical cycles were acquired.

## Later versions and predecessors

The [Xu 2021 detailed protocol](https://doi.org/10.21769/BioProtoc.4239) is an
updated method. Its original oligo section explicitly describes an inserted T in
R2R and complementary A in R2 RNA to reduce adapter dimers, giving a 35/36-nt
starter. It also offers an optional 6N UMI and optional i5 index. Those choices
are not used on this 2016 page. Its explicit E2610 adenylation or purchase of
preadenylated adapter is later evidence, not identification of the 2016 route.

The 2016 main article contrasts its ligation route with the earlier
Katibah/Shen/Zheng small-RNA/CircLigase methods. Their combined starter and cDNA
circularization are not inherited by this graph.

## Source review

The original article, EPS Figure S1 and its five oligos are recovered through
PMC Cloud. Figure S1 was visually checked and its extracted source text supports
transcription regressions. Figures S2–S9 show extraction controls and analysis
results; the data workbook contains five RNA-class count sheets. Their bounded
source checks and pending independent checks are tracked in the maintainer ledger.

The publisher-linked supplementary legends/Tables DOCX is still a concrete
retrieval obligation: its publisher route returned 403. The current NEB product and adenylated-linker FAQ have readable scientific
bodies. The linked ligation protocols currently expose title/footer shells; a
direct installed-browser attempt reached a challenge page. Their scientific
bodies remain pending rather than being treated as read from indexed extracts.
The review remains in progress, with named direct links retained and no accepted
missing-information claim.

## Bounded curation round 4 (2026-10-16)

Inspected: Published PMC4691826 JATS185–200, RNA treatments and Construction of RNA-seq libraries.

Specified 60°C/15min TGIRT extension and 12-cycle Phusion PCR. Confirmed original blocked RNA starter, preadenylated R1R ligation and dataset-specific paired75/100/150 reads; later dual-index chemistry not substituted.

Deferred: Original SupplementalFigS1 modified oligos and ordered-versus-activated R1R check, individual dataset run/index selection and two MinElute/AMPure recovery states. This pass does not establish a complete source review or publication-wide absence.
