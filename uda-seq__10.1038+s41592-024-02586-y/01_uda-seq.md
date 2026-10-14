# UDA-seq — research notes

## Sources read

- 🟢 Wang et al. 2025, DOI `10.1038/s41592-024-02586-y`, PMC article, Methods and supplements.

## Protocol boundary

🟢 UDA-seq is an indexing wrapper rather than one analyte-specific chemistry. Fixed cells or nuclei are intentionally overloaded into standard droplets and receive a droplet barcode. They remain intact when the emulsion is broken, are mixed and redistributed into wells, and receive a second well index by PCR. The barcode pair resolves cells that shared either compartment. The authors demonstrate RNA–ATAC, RNA–VDJ and RNA–CRISPR configurations; the page shows representative RNA and ATAC final libraries without claiming that one construct covers every configuration.

## Source-level oligo audit

🟢 Supplementary Table 7 (`MOESM8`, sheet “primer summary” and two 96-index sheets) prints the UDA-seq oligo families. The i5-family template is `AATGATACGGCGACCACCGAGATCT-[8-bp sample index]-ACACTCTTTCCCTACACGACGCTCTTCCGATCT`; the first row uses `GGTTTACT` and records the 3′ terminal linkage as `C-s-T`. The i7-family template is `CAAGCAGAAGACGGCATACGAGAT-[8-bp sample index]-GTCTCGTGGGCTCGG`; its first row uses `CGCCGTTT`. The same table prints `Tn5-top_ME` as `[5Phos]CTGTCTCTTATACACATCT[3ddC]`, `Tn5-bottom_Read2N` as `GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG`, a partial TSO `AAGCAGTGGTATCAACGCAGAG`, and its biotinylated counterpart `Bio-AAGCAGTGGTATCAACGCAGAG`. The Methods identify the Tn5 pair as the RNA-route i7-only transposome cassette; its modification must not be discarded.

🔴 The paper uses 10x droplet barcoding reagents and does not disclose the complete bead oligo. The exact route from all published UDA oligos into the representative Multiome library branches still needs modeling; the published table should be used before inferring any assay-dependent terminal primer.
