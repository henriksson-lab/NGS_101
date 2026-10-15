# GRO-cap

**Defining paper:** Core et al., *Nature Genetics* (2014), doi:10.1038/ng.3142. The article, methods and supplements were read from PMCID PMC4254663.

## Molecular boundary

🟢 GRO-cap begins with BrUTP-labelled global nuclear run-on RNA. It then enriches capped nascent 5′ ends: uncapped 5′ phosphates are removed, tobacco acid pyrophosphatase exposes a ligatable monophosphate only at capped ends, and a no-TAP library controls background.

🟢 The paper delegates GRO-cap library preparation to Kruesi et al. (2013). The informative coordinate is the first transcribed nucleotide, not the run-on-labelled 3′ end. Its separate PRO-seq method names TruSeq Small RNA adapters; that statement does not establish the adapter used for GRO-cap.

🟢 Kruesi et al., *eLife* 2013, doi:10.7554/eLife.00808, original PDF Methods pp. 22–23 (web-extracted lines 932–965), explicitly identifies the GRO-cap adapters as TruSeq RNA RA3 (Illumina part 15013207) and RA5 (part 15013205). After the first BrU affinity capture, RA3 is ligated to the RNA 3′ end; a second BrU capture follows. Terminator 5′ phosphate-dependent exonuclease degrades uncapped 5′ monophosphate RNA, and Antarctic phosphatase blocks remaining uncapped 5′ ends. Only then is the material split for TAP+ and no-TAP control; RA5 is ligated to the newly decapped 5′ end in the TAP+ branch, followed by a third BrU capture. The earlier graph reversed RA3/RA5 order and omitted the Terminator step; that was a repository model error.

🟢 The [official Illumina TruSeq Small RNA sequences](https://support-docs.illumina.com/SHARE/AdapterSequences/Content/TruSeq-SmallRNA.htm) identify RA3 as 5′-`TGGAATTCTCGGGTGCCAAGG`-3′ (the vendor prints T) and RA5 as 5′-`GUUCAGAGUUCUACAGUCCGACGAUC`-3′ (vendor RNA spelling). The same page prints RTP, RP1 and an RPI six-base-index family. These are kit reference sequences, not a record of the particular 2013 kit revision or chosen index.

🟡 The kit reference lists RTP 5′-`GCCTTGGCACCCGAGAATTCCA`-3′, RP1 5′-`AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA`-3′ and RPI-family 5′-`CAAGCAGAAGACGGCATACGAGAT[6-base index]GTGACTGGAGTTCCTTGGCACCCGAGAATTCCA`-3′. Kruesi prints RP1 and the RPI-family structure in its adjacent GRO-seq method (original PDF Methods p. 22, local text lines 1157–1161); Core 2008's inherited RT/PCR method names Illumina kit primers without printing individual bases. These kit sequences are candidates for the inherited GRO-cap completion, but the precise 2014 primer/index pair is unverified.

## Evidence boundary

🟡 Kruesi names RA3/RA5, and the current vendor guide supplies their sequences, but the historical kit revision and particular six-base index used by Core et al. 2014 remain unchecked. Kruesi delegates GRO-cap reverse transcription and final amplification to Core et al. 2008. Its original Science Supporting Online Material (doi:10.1126/science.1162228, pp. 18–19, web-extracted lines 507–523) uses the Illumina Small RNA RT primer with SuperScript III at 44°C then 52°C, RNase digestion, 15 Phusion PCR cycles with Illumina kit primers, and PAGE purification. That source does not identify a particular index member for the later GRO-cap libraries. The diagram marks kit-reference adapter bases as inferred for the historical experiment and leaves the selected index and final primer geometry unresolved.

## Oligo source audit

🟢 Checked the open author manuscript (`ng.3142_PMC4254663.pdf.txt`, lines 625–634) and publisher supplementary PDF (`ng.3142_41588_2014_BFng3142_MOESM77_ESM.pdf.txt`). The Methods assign TruSeq Small RNA adapters to PRO-seq, after separately assigning GRO-cap to Kruesi et al. (2013). The supplement contains figures and analyses, not a GRO-cap adapter/primer table. Kruesi's original Methods resolve the adapter names, ligation order and cap-selection chemistry. The cited Core et al. 2008 main article (doi:10.1126/science.1162228) is locally available; its original [Science Supporting Online Material](https://www.cecerelab.com/_files/ugd/594009_9144e1ea58ad4573871ec024827290f3.pdf) was independently read from a full PDF mirror because the PMC supplementary endpoint is gated. The former note incorrectly assigned the PRO-seq kit to GRO-cap; this was a repository reading error, not a published discrepancy.
