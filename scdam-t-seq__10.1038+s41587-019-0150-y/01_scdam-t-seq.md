# scDam&T-seq — research notes

## Sources read

- 🟢 Rooijers et al. 2019, DOI `10.1038/s41587-019-0150-y`, PMC full text and supplements.
- 🟢 Markodimitraki et al. 2020, DOI `10.1038/s41596-020-0314-8`, detailed protocol and supplementary oligo table.
- 🟢 Hashimshony et al. 2016 CEL-Seq2, DOI `10.1186/s13059-016-0938-8`,
  original article (predecessor of the RNA branch).

## Construct path

🟢 A Dam fusion deposits N6-methyladenine at nearby GATC sites in living cells. Cells are
sorted into wells containing a CEL-seq2 primer; its printed architecture is T7 promoter,
RA5, alternating 3-nt UMI / 4-nt cell-barcode parts, dT25 and anchoring V. Reverse
transcription is followed by second-strand synthesis in the same well.

🟢 DpnI cuts fully methylated GATC sites to blunt GA/TC half-sites. The detailed protocol
prints both strands of each DamID adapter. They contain a deliberately non-complementary
six-base fork, followed by a T7/RA5-bearing paired core, alternating UMI/barcode parts and
a terminal GA/TC blunt end. The bottom oligo is 5′ phosphorylated. Ligation of adapter
`GA` to genomic `TC` restores `GATC`; this exact junction and the unpaired fork are both
retained in the model. Supplementary Table 2 contains the complete 384-member adapter
family; the page uses the protocol's printed `TGCA`/`TATG` example.

🟢 The 2019 Rooijers Methods, “DamID adapters,” and its Supplementary Table 2
distinguish a four-base-fork adapter from the six-base-fork Dam-RING1B mESC variant.
The drawing selects the later detailed protocol's six-base-fork example; it should
not be read as the only 2019 adapter. The original CEL-Seq2 article describes a
six-base cell barcode in its own assay, whereas scDam&T-seq explicitly prints a
split eight-base barcode for its adapted primer.

🟢 CEL-seq2 and DamID use non-overlapping cell-barcode sets, but both precursors carry a
T7 promoter. Wells can therefore be pooled before joint linear amplification. The aRNA is
fragmented and reverse-transcribed using the exact
`GCCTTGGCACCCGAGAATTCCA-N6` randomhexRT primer. The printed RNA-PCR1 and indexed reverse
primer then install P5, the six-base library index and P7; their annealing suffixes and
resulting amplicon boundaries are computed on the page.
