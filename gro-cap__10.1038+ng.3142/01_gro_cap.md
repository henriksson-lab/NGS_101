# GRO-cap

**Defining paper:** Core et al., *Nature Genetics* (2014), doi:10.1038/ng.3142. The article, methods and supplements were read from PMCID PMC4254663.

## Molecular boundary

🟢 GRO-cap begins with BrUTP-labelled global nuclear run-on RNA. It then enriches capped nascent 5′ ends: uncapped 5′ phosphates are removed, tobacco acid pyrophosphatase exposes a ligatable monophosphate only at capped ends, and a no-TAP library controls background.

🟢 The paper delegates GRO-cap library preparation to Kruesi et al. (2013). The informative coordinate is the first transcribed nucleotide, not the run-on-labelled 3′ end. Its separate PRO-seq method names TruSeq Small RNA adapters; that statement does not establish the adapter used for GRO-cap.

🟢 Kruesi et al., *eLife* 2013, doi:10.7554/eLife.00808, original PDF Methods pp. 22–23 (local extracted lines 1163–1203), explicitly identifies the GRO-cap adapters as TruSeq RNA RA3 (Illumina part 15013207) and RA5 (part 15013205). After the first BrU affinity capture, RA3 is ligated to the RNA 3′ end; a second BrU capture follows. Terminator 5′ phosphate-dependent exonuclease degrades uncapped 5′ monophosphate RNA, and Antarctic phosphatase blocks remaining uncapped 5′ ends. Only then is the material split for TAP+ and no-TAP control; RA5 is ligated to the newly decapped 5′ end in the TAP+ branch, followed by a third BrU capture. The earlier graph reversed RA3/RA5 order and omitted the Terminator step; that was a repository model error.

🟢 The [official Illumina TruSeq Small RNA sequences](https://support-docs.illumina.com/SHARE/AdapterSequences/Content/TruSeq-SmallRNA.htm) identify RA3 as 5′-`TGGAATTCTCGGGTGCCAAGG`-3′ (the vendor prints T) and RA5 as 5′-`GUUCAGAGUUCUACAGUCCGACGAUC`-3′ (vendor RNA spelling). The same page prints RTP, RP1 and an RPI six-base-index family. These are kit reference sequences, not a record of the particular 2013 kit revision or chosen index.

🟢 The [Illumina TruSeq Small RNA guide, part 15004197 Rev C (2011)](https://www.bu.edu/iscf/files/2011/05/TruSeq_SmallRNA_SamplePrep_Guide_15004197_C.pdf), an original vendor guide mirrored by Boston University, lists RA3 part `15013207` and RA5 part `15013205` (local extracted lines 451–455). These match Kruesi's named part numbers. Its Table 7 (local lines 1643–1760) lists RPI1–48 index six-mers; it does not give the full RA3/RA5 strands or the index assignments for Core et al. 2014.

🟡 The kit reference lists RTP 5′-`GCCTTGGCACCCGAGAATTCCA`-3′, RP1 5′-`AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA`-3′ and RPI-family 5′-`CAAGCAGAAGACGGCATACGAGAT[6-base index]GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA`-3′. Kruesi prints RP1 and the RPI-family structure in its adjacent GRO-seq method (original PDF Methods p. 22, local text lines 1157–1161); Core 2008's inherited RT/PCR method names Illumina kit primers without printing individual bases. These kit sequences are candidates for the inherited GRO-cap completion, but the precise 2014 primer/index pair is unverified.

## Evidence boundary

🟡 The 2011 vendor guide verifies the historical component identities, and the current vendor page supplies their sequences, but complete historical ordered-strand details and Core2014 GRO-cap primer/index identity and length remain unassigned. The six-base RPI family belongs to adjacent GRO-seq and is not an established GRO-cap assignment. Kruesi delegates GRO-cap reverse transcription and final amplification to Core et al. 2008. Its original Science Supporting Online Material (doi:10.1126/science.1162228, pp. 18–19, local extracted lines 699–710) uses the Illumina Small RNA RT primer with SuperScript III at 44°C then 52°C, RNase digestion, 15 Phusion PCR cycles with Illumina kit primers, and PAGE purification. That source does not identify a particular index member for the later GRO-cap libraries. The diagram marks kit-reference adapter bases as inferred for the historical experiment and leaves the selected index and final primer geometry unresolved.

## Oligo source audit

🟢 Checked the open author manuscript (`ng.3142_PMC4254663.pdf.txt`, lines 625–634) and publisher supplementary PDF (`ng.3142_41588_2014_BFng3142_MOESM77_ESM.pdf.txt`). The Methods assign TruSeq Small RNA adapters to PRO-seq, after separately assigning GRO-cap to Kruesi et al. (2013). The supplement contains figures and analyses, not a GRO-cap adapter/primer table. Kruesi's original Methods resolve the adapter names, ligation order and cap-selection chemistry. The cited Core et al. 2008 main article (doi:10.1126/science.1162228) is locally available; its original [Science Supporting Online Material](https://www.cecerelab.com/_files/ugd/594009_9144e1ea58ad4573871ec024827290f3.pdf) was independently read from a full PDF mirror because the PMC supplementary endpoint is gated. The former note incorrectly assigned the PRO-seq kit to GRO-cap; this was a repository reading error, not a published discrepancy.

## Independent decoupling check

🟢 Rereading the original Core2014 Methods (local JATS534–549), Kruesi2013
Methods pp21–23 (local text1080–1211) and Core2008 original Science SOM
pp18–19 (local text683–710) confirms the local TAP+ order. The adapter ligations
are RA3 then RA5, with three BrU affinity captures and Terminator/phosphatase
selection before TAP. The inherited completion is RT, RNase digestion, PCR
and PAGE. The adjacent GRO-seq poly(A)/circularization and RP1/RPI declaration
and Core2014's separate PRO-seq adapter statement do not assign GRO-cap
primer strands, index length or physical library geometry.

🟢 Positive historical base evidence: Kruesi's joint “Mapping of GRO-seq and
GRO-cap reads” section (local text1205–1211) explicitly prints the cutadapt
target `TGGAATTCTCGGGTGCCAAGG`. This matches the current RA3 reference. It
supplies a historical adapter read-through target without establishing all
ordered-strand terminal modifications or the later sample-specific PCR/index
assignment. The model's conditional historical adapter representation should
be read with that positive evidence, rather than as a claim that no historical
adapter bases were printed.

The index question concerns an unassigned GRO-cap primer/index identity and
length; it must not presuppose the neighboring GRO-seq six-base RPI family.
These are model assignments to investigate, not confirmed source omissions.

## Linked study metadata

🟢 The original GEO records [GSE60456](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE60456)
and [GSE43087](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43087)
contain sample-level GRO-cap preparation descriptions, rather than just links
to the papers. The four Core2014 samples and eight Kruesi2013 samples
corroborate the RA3-first, cap-selection, RA5-second order. The precursor
metadata also prints the same RA3 trimming target as the original paper.

🟢 Original experiment records SRX682016–SRX682019 and
SRX212227–SRX212234 explicitly declare a single-read library layout on
Illumina HiSeq2000. This is a deposited read-layout statement; it does not
select a physical PCR/index primer or establish its binding geometry.
The original publisher supplement's Table1 (p29) labels its 30bp values
as **length of mapped reads**. GEO likewise labels 30bp as processing
trim, so neither value is an acquired sequencing-cycle recipe.

The paper identifies the Cornell browser-hub link as genomic tracks and
TSS results; its original URL now returns an error. The source review
records that unavailable original separately from the primer/index model
questions. No source-absence conclusion is inferred from that error.

The article separately cites Kwak for PRO-seq and Wang for GRO-seq.
Kruesi's Ingolia poly(A)/circularization reference belongs to its adjacent
GRO-seq method; its GRO-cap completion explicitly points to Core2008.
Those neighboring recipes are not assigned to this declaration.
