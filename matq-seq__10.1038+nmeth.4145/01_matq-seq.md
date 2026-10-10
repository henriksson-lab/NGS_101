# MATQ-seq — research notes

> 🟢 source transcription · 🟡 derived · 🔴 unavailable

Defining source: Sheng K et al. “Effective detection of variation in single-cell
transcriptomes using MATQ-seq.” *Nature Methods* 14, 267–270 (2017).
[doi:10.1038/nmeth.4145](https://doi.org/10.1038/nmeth.4145)

Sources read from `$CHEM_DATA/sources/matq-seq__10.1038+nmeth.4145/`:

- 🟢 the complete supplementary MATQ-seq protocol (`...MOESM227...pdf.txt`);
- 🟢 supplementary figures, especially the paired-end insert/read description;
- 🟢 supplementary data tables for context, not for reaction order.

## Reaction order

🟢 Single-cell lysate contains GAT27dT, GAT27-5N3G and GAT27-5N3T. SuperScript III
runs ten stepped-temperature cycles (8–50 °C), allowing repeated random and poly(A)
annealing. T4 DNA polymerase digests unused primers, then RNase H and RNase I remove RNA.

🟢 TdT adds dC to first-strand cDNA. `GAT21-6N3G` uses terminal GGG against that tail,
followed by ten Deep Vent exo− annealing/extension cycles. A single GAT27 primer then
amplifies for 24 cycles. `3NGAT24` converts the amplified material to more completely
double-stranded product before shearing.

🟢 All six protocol oligos are transcribed verbatim into `tools/matq_seq.py` and rendered
as raw selectable ordering text.

## Final library boundary

🟢 Supplementary Figure 21 describes paired-end 85-base reads and approximately 200-base
sheared library fragments. 🔴 The accessible source material does not print the final
Illumina adapter-kit oligos. The page therefore uses an explicitly inferred canonical
TruSeq terminal geometry and does not present it as a source transcription.
