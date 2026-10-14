# SCITO-seq

**Defining paper:** Hwang et al., *Nature Methods* (2021), doi:10.1038/s41592-021-01222-3. The open author manuscript and methods were read from PMCID PMC8643207.

## Molecular boundary

🟢 Each antibody carries a reusable 20-nt handle. A hybridized splint supplies its complement, a TruSeq Read 2 handle, antibody identity, a 10-nt staining-pool barcode and the 10x feature-capture sequence.

🟢 Cells are stained in separate pools and combined before high-loading droplet capture. Gel beads add a droplet barcode and UMI; antibody-plus-pool identity and droplet identity resolve cells combinatorially.

## Evidence boundary

🟢 The source publishes multiple splint *architectures*, including
commercial-antibody and RNA-comodality variants. For the original design it specifies,
in order, a 20-nt antibody-handle complement, 34-nt TruSeq Read 2 region, 10-nt
compound antibody/pool barcode and 22-nt 10x feature-capture region.

🟡 We rechecked the defining article and Methods, the open author manuscript,
Supplementary Tables 1–4, the bioRxiv supplement and the authors' linked SCITO-seq
analysis repository. They do not print a universal set of antibody-handle, antibody-
barcode or pool-barcode bases. Those fields are deliberately antibody/pool specific;
the final capture/indexing chemistry is inherited from the commercial 10x Chromium 3′
v3 Feature Barcoding workflow. The page therefore keeps those regions typed by role
instead of inventing a representative sequence.
