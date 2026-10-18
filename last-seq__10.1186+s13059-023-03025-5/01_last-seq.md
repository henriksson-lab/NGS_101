# LAST-seq — research notes

## Current native graph scope — 2026-10-18 batch11

Three originals were inspected: published PMC10413806 Methods150–195; original
MOESM2 XLSX TablesS2/S3; original MOESM1 DOCX FigS3/S8 legends and the rendered
FigS8 original image. Mandatory scraper:107hits in1file. Source investigation is
frozen; no further originals were opened.

🟢 Source modules form a covalent hairpin/rU-dT primer by E.coli DNA ligase,
then ureaPAGE/elution/ethanol recovery. Its mixedDNA/rU19 tail includes a physical
terminalddC residue. Cell lysis introduces independent60pgRNA carrier and ERCC
spike molecules. After evaporation, RNaseH nicks originalRNA in the RNA/dT hybrid;
Klenow extends that RNA3′end with DNA copied from the source captureprimer.
This installs a duplexDNA T7 promoter on a mixedRNA–DNA template without prior
reverse transcription. Ten independent well mixtures are pooled for the selected
route (80-well is separate), bead-recovered, directly transcribed, bead-recovered,
then reverse-transcribed with the source RandomRT primer. Originaltemplate,
physicalcaptureprimer, aRNA copies, carrier/spikes and possible genomicDNA/offtarget
copies remain independent populations or explicit unresolved residuals.

🟡 Actual random-RT ssDNA enters the first source8-cyclePCR, followed by0.9×
beads, furtherPCR and finalbeads/11µlTE. The source phase conflict is bounded:
Methods namesP5/P7 for firstPCR; TableS3 prints shortP5/P7 and NEBNextUniversal/
Index; renderedFigS8 shows completedNEBNext ends. The short primers lack full
known RT-product landing sites. The graph preserves named source reagents/order
with unavailable first-cycle geometry; it no longer silently assigns NEBNext to
firstPCR or shortP5/P7 to secondPCR. Completed endpointreference is not an incoming
PCR template. SourceUMI/cellbarcode features persist through native RNA/DNA copies.

The17bp paired linker plus11bp closing arm supplies the sourceFigS8 28bp linker
reference; the additional unpairedlinkerT belongs with the adjacentT/A context.
Exact RNaseH nick/extension positions remain unassigned, rather than fixing all
molecules to an illustrative tail length. SourceBC5/BC7 bothprintCAGGAA; the public
collision notice remains and BC2 is used as an independent printed example.
Read1/Read2=50/25; native6nt NEB index family does not establish acquired indexcycles.
Computational low-quality/emptywell filters are not physical discard.

Deferred: exact nick/extension/first-cycle geometry and PCR-stage primer assignment,
second-PCR cycles, source barcode collision and native index acquisition. No
missing-information acceptance or review-status upgrade.


> 🟢 source transcription · 🟡 derived from the listed oligos · 🔴 unavailable

Defining source: Lyu J, Chen C. “LAST-seq: single-cell RNA sequencing by direct
amplification of single-stranded RNA without prior reverse transcription and second-strand
synthesis.” *Genome Biology* 24, 184 (2023).
[doi:10.1186/s13059-023-03025-5](https://doi.org/10.1186/s13059-023-03025-5)

Sources read from `$CHEM_DATA/sources/last-seq__10.1186+s13059-023-03025-5/`:

- 🟢 PMC10413806 JATS (`.xml.txt`): Methods "LAST-seq assay" lines 147–167 and
  "Sequencing" line 191 — reaction order, conditions, identifier lengths, sequencing
  (re-read 2026-10-17).
- 🟢 Additional file 2 (`MOESM2_ESM.xlsx.txt`) Tables S2 and S3 (lines 7–134): synthetic validation oligos in S2; 80 hairpin barcode entries, dTrU tail, RandomRT and P5/P7 library primers in S3, including modifications.
- 🟢 Supplementary protocol and figures: construction of the hairpin/rU-dT primer.

## Molecular mechanism

🟢 The capture primer is assembled by E. coli DNA ligase from a self-folding hairpin
module and a 5′-phosphorylated rU-dT module ending in 3′ ddC. Its hairpin carries a PCR
handle, T7 promoter, 8-nt UMI and one of eighty 6-nt cellular barcodes. The rU tract
anneals to mRNA poly(A).

🟢 RNase H nicks the short RNA/DNA hybrid and Klenow fragment exo− extends the nick,
placing a short duplex T7 promoter beside the otherwise single-stranded original RNA.
T7 RNA polymerase then transcribes directly along that original ssRNA, producing hundreds
of antisense-RNA copies without initial RT or second-strand synthesis.

🟢 (JATS Methods "Linear amplification…", `.xml.txt` line 163) Lysates are evaporated
and treated with RNase H and Klenow exo− in NEBuffer2 (37 °C, 30 min); 10 or 80 lysates
are pooled and purified with 0.5X beads; the beads are rehydrated in T7 IVT mix and
incubated at 37 °C overnight (up to 16 h); aRNA is purified with 0.8X beads. Lysis
(line 159) is in 2 µl buffer containing 2 nM LAST-seq primer, 65 °C for 5 min. Primer
assembly (line 155): E. coli DNA ligase, 16 °C, 30 min, then 10% urea PAGE purification.

🟢 (line 167) The pooled aRNA is random-primed with `TACACGACGCTCTTCCGATCTNNNNNN`
(70 °C 2 min), reverse-transcribed with SuperScript IV (23 °C 10 min, 50 °C 15 min,
80 °C 10 min), PCR-amplified for 8 cycles with Q5 and "P5/P7 primers", purified with
0.9X beads, "further PCR amplified" and bead-purified (11 µl 0.1× TE).
Sequencing (line 191): NextSeq 550, High Output v2.5, Read 1 = 50 and Read 2 = 25;
index-read cycles are not stated.

🟡 Source ambiguity, shown on the page: Table S3 prints both short P5 (29 nt) and P7
(24 nt) "library enrichment" primers and NEBNext Universal / NEBNext1–9, 12 index
primers (3′ phosphorothioate, printed `C*T`). Computed in the model, the short primers'
3′ ends match only {{= last_seq.SHORT_SITE_NT }} nt of the RandomRT handle and hairpin PCR
handle, while the NEBNext primers have {{= last_seq.NEB_SITE_NT }}-nt sites. The historical model
assigned the 8-cycle PCR to NEBNext and the unnamed second PCR to P5/P7. The current
native graph withdraws those stage assignments: the source named P5/P7 remain an
unresolved first-cycle interpretation and the further-PCR family is unassigned. The
NEBNext endpoint family has a6nt native index; acquired index cycles remain unassigned.

🟡 Read layout derived from the hairpin sequence: Read 2 (25 cycles) reads UMI (8), cell
barcode (6) and the 11-nt hairpin closing arm `TCGTATTAGAG`; the rU-dT linker and
poly(A) follow beyond cycle 25. The final construct now includes the linker and
closing-arm complements and joins P5 to the Read 1 arm without the duplicated `ACAC`
that the earlier P5 + full TruSeq Read 1 concatenation introduced (model correction).

## Oligos retained in the model

🟢 `dTrU tail` is `/5Phos/CCACCTTTCATTCACCCTTTTTTTrU×19/3ddC/`.
The 80 barcode-specific hairpin entries printed in Supplementary Table S3 differ at their
6-nt cellular-barcode window; the model represents that window structurally rather than
choosing one cell identity. The table prints the same six-base window `CAGGAA` for
both BC5 and BC7, giving 79 distinct sequences across its 80 named entries. This is
a source-table collision, not a model-derived barcode whitelist correction.

🟡 The public final construct is assembled from the printed RandomRT, P5/P7 and hairpin
regions. It shows the barcode/UMI end as the Read-2-facing end and preserves the stated
lengths as structured features.

## Open questions

- Which primers and how many cycles were used for the second ("further") PCR.
- Whether "P5/P7 primers" in the Methods means the P5/P7-tailed NEBNext primers.
- Index-read cycles and the intended BC5/BC7 distinct barcode.
