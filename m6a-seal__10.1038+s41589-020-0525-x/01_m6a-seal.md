# m6A-SEAL — limited published-source preparation notes

Wang Y, Xiao Y, Dong S, Yu Q, Jia G. *Antibody-free enzyme-assisted chemical approach for detection of N6-methyladenosine*. Nature Chemical Biology 16, 896–903 (2020). [Published article](https://doi.org/10.1038/s41589-020-0525-x).

## Source scope

The published main HTML is a subscription preview: identity and abstract are accessible, but full Methods were not acquired. The actual published supplement MOESM1 was read at Table1, Table2 heading/selected qPCR rows, Figures1–19 captions and Note3. This is positive preparative evidence, not a complete protocol review. No full-library graph or public schematic is claimed. Published material takes precedence over automatically downloaded preprint 10.1101/790683; the preprint was not read or used. MOESM2 and data workbooks remain uninspected. The mandatory oligo scraper was run against the source bundle, and relevant source pointers were checked; hits do not establish an operational sequencing adapter.

## Distinct preparation and positive findings

🟢 Supplement Figure1 caption (lines388–391) describes sequential m6A→hm6A→dm6A conversion, followed by MTSEA attachment to the free sulfhydryl group to install a tag such as biotin. This is an enzyme-assisted chemical-labeling strategy, distinct from antibody capture or direct editing methods. It does not establish complete conversion or absolute specificity.

🟢 Figure7 caption (417–420) uses 0.5nmol FTO with 1nmol model RNA, pH7, 37°C, 5min and reports 60% hm6A formation. Figure8 (422–431) screens DTT conditions on model RNA: 200mM DTT, with pH4/3h as selected fixed conditions during the individual optimization axes. These are model-RNA experiments, not automatically the final biological-library recipe.

🟢 Figure9 (433–437) varies FTO concentration for 300ng human poly(A)+ RNA in 300µl, 5min at37°C. Figure15 (461–466) varies oxidation time using 0.2µM FTO on HEK293T poly(A)+ RNA at37°C. Do not combine these experimental scales into one inferred final reaction.

🟢 Figure12 (446–450) reports LC–MS/MS after MTSEA-biotin labeling. Figure13 (452–455) describes removal of endogenous hm6A following95°C RNA fragmentation. Figure17 (471–474) compares 5µg fragmented poly(A)+ RNA and500ng input experiments. These selected captions do not specify the complete capture/elution/library operations.

🟢 Supplement Note3 (583–617) describes independently synthesized spike-in RNAs: HiScribe T7 High Yield RNA Synthesis Kit, 500ng DNA template, 1µl T7 RNA Polymerase Mix in20µl, 37°C16h, phenol/chloroform purification; 300pg modified or unmodified spike-in is added to3µg fragmented sample before library preparation. This is a source experiment's control preparation and allocation, not a claim that every human/rice arm used that same input mass. The literal subsection calls one template “single m6A” and separately labels a capped reagent; exact transcription starts/cap installation are not modeled from prose alone.

## Oligo roles and sequencing

🟢 Supplement Table2 (156 onward) identifies qPCR enrichment primers. Table3 SELECT up/down oligos include source5phos annotations and target/control sites. These analytical primers are not assigned as sequencing-library RT/PCR adapters. Note3 supplies spike-in DNA template sequences, specified NTP pools and cap notation, with an underlined T7 promoter in the original. Their full molecules, promoter boundary and modification distributions remain unassembled; no guessed adapter sequence is adopted.

🟢 Table1 (43–44,120) names human HEK293T PE150 on HiSeq-Xten and rice SR75 on NextSeq500. Read cycles are supported by this table; index cycles, index ordering, exact sequencing primers and whole-library strand geometry remain unassigned. Human/rice applications and optimization variants are one method here.

## Deferred work

The published full article/Methods is an actual outstanding original, requested in geturl.txt. After obtaining it, resolve the operational biological-sample FTO/DTT/MTSEA recipe, physical capture/elution, library kit/version, oligo roles, stages, side products and full read geometry. Note3 directly identifies HiScribe and EpiMark N6-Methyladenosine Enrichment Kit originals for the control reagents; these are ordinary unreviewed linked sources, not established access failures. Remaining supplement originals have not been read. No missing-information claim is accepted, and related method comparisons do not establish ancestry.
