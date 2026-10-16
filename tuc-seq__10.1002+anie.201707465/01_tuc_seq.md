# TUC-seq

**Defining paper:** Riml et al., *Angewandte Chemie International Edition* (2017), doi:10.1002/anie.201707465.

## Molecular boundary

🟢 RNA synthesized during a 4sU pulse contains thiouridine. Osmium-tetroxide-mediated chemistry converts 4sU into a cytidine derivative before reverse transcription.

🟢 Newly synthesized molecules are recognized by U-to-C sequence changes distributed through the read. The conversion mark is molecular content, not an indexing sequence.

## Evidence boundary

🟢 The original publisher abstract reports detection of pulse-labelled mRNA from seven selected genes, alongside endogenous 4sU in tRNA. It does not establish a transcriptome-wide library or a universal adapter route.

🟡 Library preparation after conversion is unresolved for the 2017 experiment. The graph marks this as a warning step, and the terminal construct shows only unresolved ends around target-derived cDNA. The exact original primers, adapters, indices, and read layout cannot be inferred from that abstract.

🟢 A later same-author TUC-seq DUAL supplement describes random-hexamer reverse transcription, gene-specific PCR primers carrying barcode overhangs, and subsequent library preparation of pooled amplicons by Eurofins (section 1.12, lines 440–447 of the extracted supplement). This establishes the *later* targeted implementation, not the precise 2017 library assembly.

## Oligo source audit

🟢 The [publisher record](https://onlinelibrary.wiley.com/doi/10.1002/anie.201707465) establishes osmium-mediated 4sU-to-C conversion and targeted detection of seven mammalian mRNAs. It lists the 2.3 MB International Edition supporting-information PDF. 🔴 The full original article and both International and German supporting-information PDFs remain inaccessible through the source fetcher (publisher HTTP 403 for supplements). The [later same-author detailed chapter](https://link.springer.com/protocol/10.1007/978-1-4939-9822-7_10) exposes its abstract but gates the procedure; its publisher record also links a 2024 correction that has not been inspected. None of these absences is accepted as proof that oligos were unpublished.

🟢 `scrape_primers.py` was run on the available original-source bundle. It finds, among other hits, the Alexa647-labelled tRNA primer and gene-specific primer table in the *later* TUC-seq DUAL supplement. Those oligos cannot be assigned to the original 2017 experiment without its methods. 🔴 The 2017 oligo inventory and library-end geometry remain blocked on the original full text and supplements.
