# ChAIR

**Defining paper:** Chai et al., *Nature Methods* (2025), doi:10.1038/s41592-025-02658-7; preprint doi:10.1101/2024.05.03.592322. The open author manuscript was read from PMCID PMC12955729.

## Molecular boundary

🟢 Crosslinked nuclei undergo in-situ restriction digestion and proximity ligation, covalently joining spatially adjacent restriction fragments. Tn5 then tags accessible chromatin while contact products remain intact.

🟢 Processed nuclei enter 10x Multiome droplets. One cell barcode connects ChAIR-RNA expression, ChAIR-ATAC accessibility fragments and ChAIR-PET paired-end contact junctions.

## Evidence boundary

🟢 The published article supersedes the cited preprint as the defining paper. The page distinguishes true proximity-ligation junctions from Tn5 fragment ends rather than treating every paired DNA read as a contact.

## Oligo source audit

🟢 The published JATS full text (`s41592-025-02658-7_PMC12955729.xml.txt`, Methods “ChAIR”) gives the **bridge linker** used for proximity ligation verbatim: forward `5′-/5Phos/CGTGATATT/iBIOdT/CACGACTCT-3′` and reverse `5′-/5Phos/GAGTCGTGAAATATCACGT-3′`. Both 5′ phosphates and the internal biotinylated dT are part of the source oligo, not optional styling. Its annealed ends should be checked by computation before a base-resolved ligation state is drawn.

🟢 The same Methods then specify the Chromium Next GEM Single Cell Multiome ATAC + Gene Expression reagent kit for tagmentation, droplet barcoding and library preparation. Those kit adapters and cell barcodes are a separate vendor-defined layer, while the bridge linker above is specific to ChAIR.

🟡 The public role-level page omits a source-known custom oligo. The published PMC supplementary tables are access gated at the article-bin endpoint; the article text itself is enough to recover the bridge linker, while exact 10x kit-end treatment requires the cited kit revision.
