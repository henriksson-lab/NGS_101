# CirSeq

**Source.** Acevedo and Andino, *Nature Protocols* 9, 1760–1769 (2014), [doi:10.1038/nprot.2014.118](https://doi.org/10.1038/nprot.2014.118), PMCID PMC4418788.

Read (2026-10-17 bounded pass): PMC4418788 author-manuscript JATS text, Overview lines 28–30,
Experimental design line 46, Reagents lines 60–138, Procedure Steps 1–56 (lines 200–565); the
author-hosted eScholarship PDF was spot-checked for the same Step 18, 25, 31, 43 and 47 values.
The scraper finds no printed oligo sequences in either document.

- 🟢 Steps 1–13: Ambion AM8740 Zn2+ fragmentation (70 °C, 7.5 min), 12.5% urea-PAGE excision of 85–100-nt RNA, overnight gel elution.
- 🟢 Step 18: T4 polynucleotide kinase and T4 RNA ligase 1 act together in one T4 RNA ligase buffer/1 mM ATP mix (37 °C, 30 min); 5′ phosphorylation and circularization are therefore one pot, drawn as two graph steps.
- 🟢 Steps 22–24: random hexamers, SuperScript III (25 °C 10 min, then 42 °C); RNase H is added after 2 min at 42 °C, then 30 min more. Rolling-circle RT creates tandem cDNA copies of one circular RNA fragment.
- 🟢 Steps 25–35: NEBNext second-strand synthesis module (16 °C, 2.5 h), NEBNext end repair (20 °C, 30 min), Klenow (3′→5′ exo−) dA-tailing (37 °C, 30 min), NEBNext Quick T4 DNA ligase with TruSeq indexed adapter (20 °C 15 min or 16 °C overnight); phenol–chloroform/ethanol cleanups between steps.
- 🟢 Steps 36–45: 10% urea-PAGE excision of 450–600-base adapter-ligated DNA; at least 435 bases = 135 adapter bases plus at least 300 bases of tandem-repeat insert.
- 🟢 Steps 46–53: Phusion PCR with the TruSeq PCR primer cocktail, cycles 2–16 (15 amplification cycles) at 98/65/72 °C, then 7.5% native PAGE purification.
- 🟢 Line 46 recommends 300-nt reads on MiSeq (MiSeq v2 300-cycle kit, line 138); HiSeq 1500/2500 Rapid mode with SR cluster kits is an alternative. 🟡 A single 300-nt read is derived from these kit and read-length statements; no paired or index read is described. Consensus uses repeated copies inside one read, without a barcode.
- 🔴 The protocol does not describe how the NEBNext module primes the second strand (the model no longer declares an unsupported "random second-strand primer"), and does not print the PCR primer-cocktail sequences; that PCR is drawn without primer geometry.
- 🟡 Open: the drawn canonical TruSeq single-index adapter segments (from `lib/illumina`) total fewer than the 135 adapter bases stated at Step 43. The kit-original adapter sequences were not inspected, so this is an unresolved model/source difference, not a published error.

The named FC-121-4001 kit is the Illumina TruSeq Nano DNA LT Set A kit, which has twelve
six-base index choices ([Illumina guide, kit contents and index table](https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_truseq/truseqnanodna/truseq-nano-dna-library-prep-guide-15041110-d.pdf),
pp. 28–29, 34). This 2015 revision postdates the 2014 paper; its revision history traces
the kit guide to May 2013. The protocol does not assign a particular Set A member or
specify index cycles for the reported run. An adapter Read 2 landing site does not
establish that a paired Read 2 was acquired.

Mapping the consensus to a known reference resolves the RNA circle's 3′→5′ ligation junction.
