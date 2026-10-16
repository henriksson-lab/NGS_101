# Start-seq

**Defining paper:** Nechaev et al., *Science* (2010), doi:10.1126/science.1181421. The open author manuscript was read from PMCID PMC3435875.

## Molecular boundary

🟢 The method isolates short capped nuclear RNAs associated with promoter-proximal Pol II. The published analysis uses the 5′ RNA end to locate transcription initiation and a purpose-built second library to read RNA 3′ ends and locate pausing.

🟢 The informative molecules are short Start-RNAs rather than fragments sampled uniformly from total RNA.

## Evidence boundary

🟢 The accessible manuscript states that its custom 3′-end adapters remain compatible with commercial Illumina sequencing primers. The original [GEO GSE18643 sample records](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18643) give reaction order for both endpoint libraries. The supplementary file that may give their custom oligo bases remains behind the PMC download gate, so the schematic shows the two endpoints without inventing bases.

## Oligo source audit

🟢 The open author manuscript (`science.1181421_PMC3435875.xml.txt`, lines 38–44) says it isolates short, nuclear, capped RNA and designs new adapters for direct 3′-end sequencing. It cites figures S2–S3 and Table S1 but does not say Table S1 is an oligo table. The [5′-end GEO sample GSM463298](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM463298) and [3′-end sample GSM463301](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM463301) are original author-deposited methods: both gel-select 25–100-nt RNA, treat with 5′ polyphosphatase then Terminator exonuclease to retain capped RNA, ligate a 3′ linker, gel-purify, treat with alkaline phosphatase then tobacco acid pyrophosphatase, and ligate a 5′ linker. The 5′ library uses standard Illumina small-RNA linkers; the 3′ library uses two custom linkers and a custom RT primer. This is the source for the two distinct graph branches and corrects an earlier model error that used the observed 25–65-nt pause interval as the size-selection window.

🟡 The official [Illumina Small RNA Sample Preparation Guide Rev D (August 2010)](https://sapac.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_legacy/SmallRNA_SamplePrep_Guide_1004239_D.pdf) identifies SRA 5′/3′ adapters and RT primer by reagent part number, but prints no bases in the fetched text. Its standard workflow ligates the 5′ adapter first; the original GEO methods explicitly use a modified 3′-then-5′ order. Rev D was issued after the original experiment and does not establish the authors' exact kit revision or custom oligo sequences.

🔴 The linked original `NIHMS401007-supplement-supp_material.pdf` is unavailable through the PMC route; Science's alternate supplement route returned 403. The later Henriques 2013 paper ([10.1016/j.molcel.2013.10.001](https://doi.org/10.1016/j.molcel.2013.10.001), JATS lines 142–148) reuses Nechaev's 3′ scRNA library procedure and does not print the custom adapter bases in its accessible article. Its own supplement is gated too. Exact custom 3′-library adapters, RT primer and any sample-index identity remain blocked pending these sources; the page no longer asserts an index.
