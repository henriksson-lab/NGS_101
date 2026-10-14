# scMT-seq — source intake

Defining paper: Hu et al., *Genome Biology* (2016), DOI
[10.1186/s13059-016-0950-z](https://doi.org/10.1186/s13059-016-0950-z).

🟢 The physically separated fractions enter transcriptome and bisulfite-methylome
branches. scMT-seq remains distinct from scM&T-seq, a separate named protocol.

## Sources and inherited chemistries

- 🟢 Hu et al. explicitly state that cytoplasmic cDNA is amplified by Smart-seq2 (with SuperScript III substituted for SuperScript II), followed by Nextera tagmentation and library construction.
- 🟢 The nuclear fraction follows the previously published scRRBS workflow: MspI digestion, end repair/dA-tailing, indexed premethylated adapter ligation, bisulfite conversion and two PCR rounds.
- 🟢 Smart-seq2 oligo-dT, TSO and ISPCR bases are taken from the Picelli supplementary protocol through the shared `smartseq2` model. The scRRBS adapter and final-library structures are taken through the shared `rrbs` model.

The defining paper does not introduce a new oligo set for either branch. The earlier page's generic post-bisulfite random-primer branch was therefore incorrect for scMT-seq and has been replaced by the inherited scRRBS chemistry.
