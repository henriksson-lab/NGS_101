# LAST-seq — research notes

> 🟢 source transcription · 🟡 derived from the listed oligos · 🔴 unavailable

Defining source: Lyu J, Chen C. “LAST-seq: single-cell RNA sequencing by direct
amplification of single-stranded RNA without prior reverse transcription and second-strand
synthesis.” *Genome Biology* 24, 184 (2023).
[doi:10.1186/s13059-023-03025-5](https://doi.org/10.1186/s13059-023-03025-5)

Sources read from `$CHEM_DATA/sources/last-seq__10.1186+s13059-023-03025-5/`:

- 🟢 PMC full text/JATS: reaction order, identifier lengths, sequencing.
- 🟢 Supplementary Table S2: oligos verbatim, including modifications.
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

🟢 The pooled aRNA is random-primed with
`TACACGACGCTCTTCCGATCTNNNNNN`, reverse-transcribed, and amplified with P5/P7 primers.
Sequencing used NextSeq 550, Read 1 = 50 and Read 2 = 25.

## Oligos retained in the model

🟢 `dTrU tail` is `/5Phos/CCACCTTTCATTCACCCTTTTTTTrU×19/3ddC/`.
The 29 barcode-specific hairpin oligos printed in Supplementary Table S2 differ at their
6-nt cellular-barcode window; the model represents that window structurally rather than
choosing one cell identity.

🟡 The public final construct is assembled from the printed RandomRT, P5/P7 and hairpin
regions. It shows the barcode/UMI end as the Read-2-facing end and preserves the stated
lengths as structured features.
