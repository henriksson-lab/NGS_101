# NTT-seq — source notes

Defining paper: Stuart et al., *Nature Biotechnology* (2023), [doi:10.1038/s41587-022-01588-5](https://doi.org/10.1038/s41587-022-01588-5).

- 🟢 PMC article, Methods and Extended Data Table 2 were read.
- 🟢 Nanobody–Tn5 fusion proteins recognize primary antibodies and transfer one of six printed 8-nt MEDSA modality barcodes; MEDSB supplies the other Tn5 end.
- 🟢 Extended Data Table 2 prints MEDSA_1 as `TCGTCGGCAGCGTC-GGATTGCT-GCGATCGAGGACGGC-AGATGTGTATAAGAGACAG` and MEDSB as `GTCTCGTGGGCTCGG-AGATGTGTATAAGAGACAG`. The diagram splits these into S5/S7, modality barcode, connector and mosaic end rather than abbreviating the transposome.
- 🟢 Each adapter is shown annealed only across its terminal 19-nt mosaic end to the 5′-phosphorylated mosaic-end reverse oligo. The transferred MEDSA_1–DNA–MEDSB intermediate is shown before droplet processing.
- 🟢 Single-cell NTT-seq then follows Chromium single-cell ATAC; the paper prints custom Read 1 and i5 primers.
- 🟢 Reported sequencing is i5 38 bp, i7 8 bp, Read 1 60 bp and Read 2 60 bp; the PBMC run without surface proteins used 50+50 genomic reads.
- 🟡 One representative MEDSA barcode is drawn. Its feature role is modality, not cell identity.
- 🟢 The linked 10x Chromium Next GEM Single Cell ATAC v1.1 guide CG000209 Rev G, Appendix A1 (lines 3729–3829 in the extracted copy), prints the gel-bead P5–N16–partial Read1N oligo and SI-PCR primer scaffolds. The cell barcode and sample index are variable windows. The representative strand-resolved PCR binding geometry has not yet been reconstructed.

The schematic represents the chromatin library branch. 🟢 The original Methods (JATS lines 110–112 and 128–136; bioRxiv v1 lines 603–634) also describe an optional TotalSeq-A surface-protein branch: an N9V/poly(dT) bridge oligo with `/3InvdT/` is added to droplet barcoding, short antibody tags are recovered from the SPRI supernatant, and a separate P5/RPI-x PCR makes the tag library. The [author's library-structure page](https://ntt-seq.com/html/nttseq_lib_struct) prints the chromatin branch's custom Read 1 and Index 2 sequencing primers and final duplex. A complete molecule model for the separate antibody-tag product remains open.
