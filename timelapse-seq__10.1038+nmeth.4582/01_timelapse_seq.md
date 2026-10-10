# TimeLapse-seq

**Defining paper:** Schofield et al., *Nature Methods* (2018), doi:10.1038/nmeth.4582. The open article and supplementary protocol were read from PMCID PMC5927882.

## Molecular boundary

🟢 Cells incorporate 4-thiouridine (4sU) into newly synthesized RNA during a pulse. Oxidation followed by 2,2,2-trifluoroethylamine treatment changes the pairing behavior of 4sU before reverse transcription.

🟢 Treated 4sU positions are recovered as T-to-C mismatches in ordinary RNA-seq reads. The mutation is the temporal label; it is not a UMI or sample barcode.

## Evidence boundary

🟡 TimeLapse chemistry is compatible with several downstream RNA-seq preparations. The schematic therefore fixes the recoding reaction and readout, but labels kit-dependent library ends by role.
