# SNARE-seq — research notes

## Sources read

- 🟢 Chen et al. 2019, DOI `10.1038/s41587-019-0290-0`, PMC full text, Supplementary Protocol and Supplementary Table 4.
- The Teichlab SNARE-seq construct page was used for orientation only; it is not source evidence.

🟢 **Published table cross-reference error:** Chen et al. Methods (original PMC
JATS line 64) directs the Nextera-R1-rc-polyA splint to “Table S3”. The actual
Table S3 (`sup_tab3.xlsx.txt`, lines 1–4) contains differential-accessibility
results; it does not list oligos. Original Table S4 (`sup_tab4.xlsx.txt`, sheet
“Oligonucleotide sequences”, lines 1–27) prints `/5Phos/`–`CTGTCTCTTATACACATCTGACGCTGCCGACG`–A31 for that splint, plus blockers,
TSO, PCR, index and custom Read 1 primers. The model uses the Table S4 strand.
There is no disagreement about the printed bases; the error is the Methods
table number. The particular index-primer member used for each library is not
established by this correction.

## Construct path

🟢 Nuclei are tagmented before droplets. Supplementary Table 4 prints
Nextera-R1-rc-polyA as `/5Phos/`–ME reverse complement–S5 reverse
complement–A31. Its complete 33-nt Nextera arm therefore anneals across an S5+ME
transposition end, while A31 pairs with the bead's terminal dT30. The page derives this
three-strand capture intermediate from those exact oligos rather than drawing a generic
"splint". 🟢 A 12-base cell barcode is therefore shared between RNA and accessibility
products; the adjacent 8-base UMI belongs to the RNA readout.

🟢 The combined reverse-transcription/ligation mix contains Hemo KlenTaq, T7 ligase and
reverse transcriptase. The captured chromatin intermediate is gap-filled, and the
splint-side nick can be sealed because Nextera-R1-rc-polyA is 5′ phosphorylated. The
other nick is not shown as covalently sealed: 🟡 the ordered S5 strand in the source
oligo table has no stated 5′ phosphate, while a completed bottom strand could support PCR.

🟢 The common amplification uses Nextera-R2 with TSO-PCR. Their exact sequences select
the asymmetric S5/S7 chromatin product and retain the bead barcode. The RNA branch uses
template switching and TSO-PCR sites at both ends. The products are then split: RNA is
TSO-enriched and Nextera-XT-finished, whereas chromatin receives P5xx-TSO/Ad2.x indexes.
Supplementary Table 4 is the authority for all modified and amplification oligos shown.

🟢 The same table also prints six `P50x-Tso` i5 primers, six `Ad2.x` i7 primers and the `Read1CustomSeqB` sequencing primer. The page now includes the first member of each index family and the exact custom Read 1 primer; `P501-Tso` retains both terminal phosphorothioate marks (`T*A*C`). Which pair of index members was selected varies by library and is not claimed by the representative examples.
