# NTT-seq — source notes

Defining paper: Stuart et al., *Nature Biotechnology* (2023), [doi:10.1038/s41587-022-01588-5](https://doi.org/10.1038/s41587-022-01588-5).

- 🟢 PMC article, Methods and Extended Data Table 2 were read.
- 🟢 Nanobody–Tn5 fusion proteins recognize primary antibodies and transfer one of six printed 8-nt MEDSA modality barcodes; MEDSB supplies the other Tn5 end.
- 🟢 Extended Data Table 2 prints MEDSA_1 as `TCGTCGGCAGCGTC-GGATTGCT-GCGATCGAGGACGGC-AGATGTGTATAAGAGACAG` and MEDSB as `GTCTCGTGGGCTCGG-AGATGTGTATAAGAGACAG`. The diagram splits these into S5/S7, modality barcode, connector and mosaic end rather than abbreviating the transposome.
- 🟢 Each adapter is shown annealed only across its terminal 19-nt mosaic end to the 5′-phosphorylated mosaic-end reverse oligo. The transferred MEDSA_1–DNA–MEDSB intermediate is shown before droplet processing.
- 🟢 Single-cell NTT-seq then follows Chromium single-cell ATAC; the paper prints custom Read 1 and i5 primers.
- 🟢 Reported sequencing is i5 38 bp, i7 8 bp, Read 1 60 bp and Read 2 60 bp; the PBMC run without surface proteins used 50+50 genomic reads.
- 🟡 One representative MEDSA barcode is drawn. Its feature role is modality, not cell identity.
- 🔴 The paper names the 10x scATAC workflow but does not disclose the complete gel-bead or SI-PCR oligos. Their sequence and PCR binding geometry are therefore unavailable; the inferred handle on the schematic exists only to preserve the published molecule order.
