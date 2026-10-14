# smol-seq

**Defining paper:** Tan and Fraser, *Nature Biotechnology* (2026), doi:10.1038/s41587-025-02554-7. The open paper was read on the publisher site.

## Molecular boundary

🟢 Each structure-switching aptamer recognizes one small-molecule target and retains a DNA barcode identifying that sensor. Target binding changes aptamer structure and releases the barcode.

🟢 Released barcodes are PCR-indexed and counted by DNA sequencing. The sequenced object identifies the aptamer and hence the analyte; the metabolite itself is never sequenced.

## Evidence boundary

🟢 The paper demonstrates multiplexing with target-specific aptamer/barcode pairs. Those pairs are sensor definitions rather than universal sequencing adapters, so the schematic encodes their feature-barcode role.

## Oligo source audit

🟢 The published `s41587-025-02554-7_supp_41587_2025_2554_MOESM2_ESM.xlsx`, sheet `SuppData1 Oligonucleotides`, gives actual ligand-binding oligo (LBO) and short release oligo (SRO) sequences, with 5′ biotin and other synthesis modifications retained. For example, the ATP LBO is `/5Biosg/TAAGATGACTCGCGAGAGCCCTGGGGGAGTATTGCGGAGGAAGGGCTCTCGC`; its barcoded ATP SRO examples carry distinct identity windows and a downstream `N5-7` random region. The supplement lists comparable sensor pairs for mefloquine and piperaquine in the barcode-sequencing experiment, as well as fluorescent or T7 reporters for assay characterization. The fluorescent `/36-FAM/` and spacer `/iSp18/` variants must not be silently substituted for the DNA-sequenced barcodes.

🟢 The same sheet lists the library primers: round 1 `Reverse P5 primer` `ACACTCTTTCCCTACACGACGCTCTTCCGATC`, three analyte-specific forward primers with the shared `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT` handle, and round 2 primers `AATGATACGGCGACCACCGAGATCTACAC [i5 index] ACACTCTTTCCCTACACGACG` and `CAAGCAGAAGACGGCATACGAGAT [i7 index] GTGACTGGAGTTCAGACGTG`. These are source-known handles; i5/i7 values are selected during library preparation.

🟡 The target-specific LBO/SRO sequences are published experimental sensor choices. An exact molecular model should select a demonstrated sensor/barcode and distinguish it from the other reporters; a single generic aptamer cannot encode all of them.
