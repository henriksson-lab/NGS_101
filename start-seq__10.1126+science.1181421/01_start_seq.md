# Start-seq

**Defining paper:** Nechaev et al., *Science* (2010), doi:10.1126/science.1181421. The open author manuscript was read from PMCID PMC3435875.

## Molecular boundary

🟢 The method isolates short capped nuclear RNAs associated with promoter-proximal Pol II. The published analysis uses the 5′ RNA end to locate transcription initiation and a purpose-built second library to read RNA 3′ ends and locate pausing.

🟢 The informative molecules are short Start-RNAs rather than fragments sampled uniformly from total RNA.

## Evidence boundary

🟡 The accessible manuscript states that the new 3′-end adapters remain compatible with commercial Illumina sequencing primers, but the supplementary oligo table is behind the PMC download gate. The schematic therefore shows the two endpoint libraries without inventing bases.

## Oligo source audit

🟢 The open author manuscript (`science.1181421_PMC3435875.xml.txt`, line 44) explicitly says that new RNA adapters for direct 3′-end sequencing are detailed in figures S2–S3 and Table S1. The original 5′-end library and the purpose-built 3′-end library therefore need separate adapter, ligation and sequencing-primer inventories. 🔴 The linked `NIHMS401007-supplement-supp_material.pdf` is unavailable through the PMC automated route (recorded as a manual row in the source manifest); its Table S1 has not been read. Neither library's exact oligos can be marked source-verified from the main manuscript alone.
