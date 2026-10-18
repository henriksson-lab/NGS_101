# Quartz-Seq2 — research-note pointer

## Current native graph scope — 2026-10-18 batch10

Two originals were inspected: published PMC5845169 JATS Methods232–258
(barcode design, sorting/ERCC, native WTA, ligation-library branch and sequencing)
and original Additionalfile5/MOESM5 TableS4 XLSX (v3.2 heading/first-member rows
389–397 and all relevant related-oligo rows near1929–1945). They remain in the
existing external Quartz-Seq source collection, not a reconstructed repository
source directory. Mandatory scraping of that collection found636hits in11files;
original pointers were verified. The dedicated-bundle oligo guard skips while
these optional source-transcription tests read the explicitly identified originals.

🟢 Selected v3.2: M–AC–15-base well barcode–N8DNA UMI–dT24 RT primer;
SuperScriptIII35°C5min/50°C50min/70°C15min →physical pooling of independent
well mixtures →DCC5 columns →TdT+RNaseH37°C75s/65°C10min →physical
45µl-to-four-approximately11µl aliquoting →M-dT23 tagging with0.2°C/s
increment40→68°C →11-cycle gM PCR in each tube →pool →MinElute →AMPure
→CovarisLE220 →DCC5 →combined KAPA end-repair/A-tail37°C60min/65°C30min
→truncated adapter ligation20°C15min →AMPure →P5-gMac/TPC2 selective
PCR8cycles →AMPure. The main Methods11-cycle value is used; no separate
v3.2 nine-cycle claim is inherited from the historical combined note.

🟡 The actual second DNA strand, plus original first strands and residual RNA,
enters WTA PCR. Native gM has no amine modification assigned from TableS4.
The source LT3′T strand covalently joins the opposite fragment5′end. That
native LT-plus-opposite-cDNA strand, including its terminal original end-prep
dA, is the selective-PCR template. The separate RT-end fragment retains its
source gM primer lineage and upstream ligated LT strand outside the internal
gM PCR boundary. Neither strand has a prematurely copied adapter arm or completed
P5-gMac extension before PCR. Short-adapter5′ligation/phosphorylation geometry is
unassigned; lack of a printed modification is not an accepted absence claim.
The two P5-gMac and one TPC2 phosphorothioate bonds remain primer-derived bonds.

The representative well is original eMDRT0001; other pooled wells retain their
own source barcode members and N8UMIs as separate molecules. Independent ERCC
RNA/copies are not concatenated with cell RNA/copies. RNA hydrolysates, original
DNA, internal/tagging-end/offtarget fragments and oligo residuals have unresolved
class-specific recovery fates; nonamplification and computational filtering are
not physical waste. The selected v3.2 acquisition is R1/I1/R2=23/6/63; v3.1
22/6/64–118 and alternative NexteraXT preparation are separate source branches.

The existing published14/15-mer heading discrepancy is retained below/on the
public page. Deferred: MOESM1 development figures/notes, buffer TableS5 and
other supplements, direct KAPA/vendor geometry, alternative NexteraXT branch,
and exact molecule-class recovery efficiencies. This capped pass does not accept
absent information, assert a full source review or upgrade historical status.


Quartz-Seq2 is documented from Sasagawa *et al.*, 2018,
doi:[10.1186/s13059-018-1407-3](https://doi.org/10.1186/s13059-018-1407-3).

The detailed source record was originally written jointly for the two related protocols
and is retained as provenance in the
[combined Quartz-Seq / Quartz-Seq2 research note](../quartz-seq__10.1186+gb-2013-14-4-r31/01_quartz-seq.md#b-quartz-seq2-2018).
The public schematic itself uses only the Quartz-Seq2 paper and its supplementary oligo
table for its molecular construct.

## Published Table S4 discrepancy

🟢 Additional file 5, Table S4 (MOESM5 XLSX text lines 389–393), heads the v3.2
sheet as “SeqLv15mer” but its next line says “14 mer.” The first printed v3.2
barcode, eMDRT0001 `ATCAATCCTATCTGC`, is 15 bases; the defining article's
Methods (“Design of cell barcodes”) likewise specifies 15 nt for v3.2.
🟡 The schematic models a 15-base v3.2 barcode because the actual oligos and
Methods agree. The 14-mer line appears to be a table heading error; a formal
correction notice was not found in the inspected sources.

## Bounded curation round 2 (2026-10-16)

Inspected: Original PMC5845169 JATS60,72,236–240,256–268, stored in Quartz-seq source bundle.

Added 11-cycle suppression PCR. Source confirms selected v3.2 read mask23/6/63; v3.1 is a different mask. Public notice qualifies RNase H activity order in the tailing mixture.

Deferred: Original primer tables and supplements; combined RNase H/TdT graph activity; dedicated source bundle absent, original read from existing external Quartz-seq collection. This capped pass does not establish a complete source review or publication-wide absence.
