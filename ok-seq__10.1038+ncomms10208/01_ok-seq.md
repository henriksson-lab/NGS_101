# OK-seq — research notes

## Sources read

- 🟢 Petryk et al. 2016, DOI `10.1038/ncomms10208`, PMC full text and supplement.
- 🟡 Wu et al. 2023, DOI `10.1038/s41596-022-00793-5`, detailed protocol cited by the earlier reconstruction; the source fetcher could not retrieve its full text for this audit.

## Construct path

🟢 EdU-labelled genomic DNA is denatured and short single strands are recovered by sucrose-gradient size fractionation. EdU is click-biotinylated, RNA primers are hydrolysed and DNA 5′ ends are phosphorylated. 🟢 Two double-stranded adapters with six-base random overhangs ligate to the single-stranded fragment; biotin capture precedes PCR. The defining human paper explicitly publishes the four adapter oligos shown on the page.

🟢 The defining Petryk paper itself prints the complete forward PCR primer, the reverse primer carrying a six-base i7 index, and the A1 Watson oligo used as the custom sequencing primer. Its forward primer begins `ATGATAC…`; the earlier model had silently prepended an extra `A`, which has been corrected. The reverse PCR primer completes the final base of the A2/TruSeq Index-1 binding site; the schematic therefore distinguishes the ligated pre-PCR adapter from the completed library rather than treating them as identical.
