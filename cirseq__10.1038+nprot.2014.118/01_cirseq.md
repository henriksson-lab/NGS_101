# CirSeq

**Source.** Acevedo and Andino, *Nature Protocols* 9, 1760–1769 (2014), [doi:10.1038/nprot.2014.118](https://doi.org/10.1038/nprot.2014.118), PMCID PMC4418788.

- 🟢 Steps 1–18: Zn2+ fragmentation, denaturing-PAGE selection of 85–100-nt RNA, 5′ phosphorylation, and T4 RNA ligase 1 circularization.
- 🟢 Steps 19–24: random-primed rolling-circle reverse transcription creates tandem cDNA copies of one circular RNA fragment.
- 🟢 Steps 25–53: second-strand synthesis, blunting, dA-tailing, TruSeq indexed-adapter ligation, size selection, and PCR.
- 🟢 Adapter-ligated molecules of 450–600 bases contain 135 adapter bases and at least 300 bases of tandem-repeat insert.
- 🟢 The described run used single-end 300-cycle MiSeq sequencing. Consensus uses repeated copies inside one read, without a barcode.

The named FC-121-4001 kit is the Illumina TruSeq Nano DNA LT Set A kit, which has twelve
six-base index choices ([Illumina guide, kit contents and index table](https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_truseq/truseqnanodna/truseq-nano-dna-library-prep-guide-15041110-d.pdf),
pp. 28–29, 34). This 2015 revision postdates the 2014 paper; its revision history traces
the kit guide to May 2013. The protocol does not assign a particular Set A member or
specify index cycles for the reported run. An adapter Read 2 landing site does not
establish that a paired Read 2 was acquired.

Mapping the consensus to a known reference resolves the RNA circle's 3′→5′ ligation junction.
