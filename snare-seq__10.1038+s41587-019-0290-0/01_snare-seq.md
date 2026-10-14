# SNARE-seq — research notes

## Sources read

- 🟢 Chen et al. 2019, DOI `10.1038/s41587-019-0290-0`, PMC full text, Supplementary Protocol and Supplementary Table 4.
- 🟢 The Teichlab SNARE-seq construct page, checked as a secondary reconstruction.

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
other nick is not shown as covalently sealed: the source reconstruction notes that its S5
end lacks a phosphate and that the completed bottom strand is sufficient for PCR.

🟢 The common amplification uses Nextera-R2 with TSO-PCR. Their exact sequences select
the asymmetric S5/S7 chromatin product and retain the bead barcode. The RNA branch uses
template switching and TSO-PCR sites at both ends. The products are then split: RNA is
TSO-enriched and Nextera-XT-finished, whereas chromatin receives P5xx-TSO/Ad2.x indexes.
Supplementary Table 4 is the authority for all modified and amplification oligos shown.

🟢 The same table also prints six `P50x-Tso` i5 primers, six `Ad2.x` i7 primers and the `Read1CustomSeqB` sequencing primer. The page now includes the first member of each index family and the exact custom Read 1 primer; `P501-Tso` retains both terminal phosphorothioate marks (`T*A*C`). Which pair of index members was selected varies by library and is not claimed by the representative examples.
