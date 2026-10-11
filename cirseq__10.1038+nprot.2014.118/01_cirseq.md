# CirSeq

**Source.** Acevedo and Andino, *Nature Protocols* 9, 1760–1769 (2014), [doi:10.1038/nprot.2014.118](https://doi.org/10.1038/nprot.2014.118), PMCID PMC4418788.

- 🟢 Steps 1–18: Zn2+ fragmentation, denaturing-PAGE selection of 85–100-nt RNA, 5′ phosphorylation, and T4 RNA ligase 1 circularization.
- 🟢 Steps 19–24: random-primed rolling-circle reverse transcription creates tandem cDNA copies of one circular RNA fragment.
- 🟢 Steps 25–53: second-strand synthesis, blunting, dA-tailing, TruSeq indexed-adapter ligation, size selection, and PCR.
- 🟢 Adapter-ligated molecules of 450–600 bases contain 135 adapter bases and at least 300 bases of tandem-repeat insert.
- 🟢 The described run used single-end 300-cycle MiSeq sequencing. Consensus uses repeated copies inside one read, without a barcode.

Mapping the consensus to a known reference resolves the RNA circle's 3′→5′ ligation junction.
