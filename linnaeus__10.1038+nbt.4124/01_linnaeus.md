# LINNAEUS — RFP scars and transcriptome readout

Spanjaard et al. (2018), *Simultaneous lineage tracing and cell-type identification using CRISPR–Cas9-induced genetic scars*, [10.1038/nbt.4124](https://doi.org/10.1038/nbt.4124).

## Sources and conditional version

Selected originals: **accepted author manuscript PMC5942543**, Methods **Scar detection in bulk samples** and **Transcriptome and scar detection in single cells** (82–90), plus source recorder overview; **published publisher SI MOESM11**, Table 1, rendered original PDF p. 41, checked for actual primer row assignment/line wraps; and final publisher HTML, whose identity/abstract are accessible but whose body is explicitly a **subscription preview**, without final Methods. The manuscript SI is separately cached and was not substituted for the publisher SI. The standard fetcher and full-bundle oligo scraper were run, then the displayed primers were checked in the rendered published Table 1.

Operational graphs are conditional accepted-manuscript routes; the published SI primer inventory is separately qualified. Matching named components does not establish complete article-version equivalence. No current Chromium recipe is substituted for its unreviewed historical original. These capped checks do not establish a complete source/linked-material or missing-information review.

## Recorder context

🟢 The source uses Zebrabow M with **16–32 independently integrated RFP reporters**, scarred by a guide targeting RFP. This differs from a single compact multi-target GESTALT array. Reporter creation/scarring occurs upstream of library preparation. The model records inherited scars as structured lineage barcodes, independent of a present-cell identifier, transcript UMI and sample index. It does not predict scar sequences, integration positions, one unique scar per cell, or indel frequencies.

## Single-cell route

🟢 Historical Chromium PN-120233 encapsulates cells and barcoded beads, followed by lysis and polyT RT. Cell barcode, UMI and adapter components are stated, but their exact lengths/order and historical kit guide remain unassembled. After pooling and initial amplification, the **whole cDNA mixture is split in half**. The first half is fragmented and continued into a conventional transcriptome library. The second half is **unfragmented**, undergoing two RFP-targeted PCRs with nested forward primers and a reverse adapter primer. The allocation is not a physical purification of RFP versus other transcripts; source reaction mixtures remain heterogeneous until an actually described separation.

🟢 The source reports NextSeq 500 paired-end 75 and HiSeq 2500 paired-end 100, without resolving this generic model's exact sample/run-primer mapping. Scar and transcript reads retain matching cell/UMI identity. Computational scar truncation to 75 nt, error filtering and lineage-tree reconstruction do not physically trim DNA or remove molecules during preparation.

## Bulk routes

🟢 Genomic DNA: heat source samples in 50 mM NaOH at 95°C for 20 min, neutralize with one-tenth volume of 1 M Tris-HCl pH 8.4, and take 20 µL for barcoded RFP PCR. Neutralization is not a purification. Pool sample-barcoded PCR products, perform AMPure cleanup, then add Illumina adapters in a second PCR. Two symbolic samples/primer members show physical pooling without assigning particular experimental animals to those oligos. Untransferred lysate remains unresolved, not discarded by assumption.

🟢 Bulk RNA: source TRIzol isolation is followed by isopropanol precipitation, two 75% ethanol washes and air drying. Source polyT/SuperScript II RT at 42°C for 2 h precedes the same specifically cited bulk RFP PCR/pooling/cleanup/adapter-PCR route. It does not acquire a single-cell barcode or UMI merely because the single-cell assay uses them. Exact TRIzol phase chemistry and polyT primer length remain unassembled.

## Published primer inventory

| Published Table 1 input | Sequence, 5′→3′ |
|---|---|
| Single-cell RFP forward 1 | `{{= chem.RFP1 }}` |
| Single-cell RFP forward 2 | `{{= chem.RFP2 }}` |
| Single-cell reverse | `{{= chem.REVERSE }}` |
| Bulk reverse | `{{= chem.BULK_REVERSE }}` |
| Bulk forward 1 | `{{= chem.bulk_forward('a').top() }}` |
| Bulk forward 2 | `{{= chem.bulk_forward('b').top() }}` |

The forward bulk primers differ in their inline sample barcode. The displayed two members carry distinct source eight-base windows; no source-specific sequencer i5/i7 orientation is inferred from those inline sample tags. Canonical matching Read1/Read2 sequence regions are imported, without assuming a canonical full library recipe. The sgRNA-synthesis primers in the same table concern recorder editing and are not sequencing-library primers.

## Boundaries and remaining work

No scGESTALT/inDrops or CARLIN recipe is imported. Related recording concepts alone do not establish exact method ancestry. Full final-version Methods, historical Chromium chemistry/version, recorder/primer/product/end geometry, bulk adapter-PCR oligos, exact run/primer geometry, remaining original/linked sources and independent review remain pending. All previews are private and incomplete; no absent-information claim is accepted.
