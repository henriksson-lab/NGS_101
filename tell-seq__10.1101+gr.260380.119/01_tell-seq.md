# TELL-seq — research notes

## Defining source

Chen Z, Pham L, Wu T-C, et al. *Ultralow-input single-tube linked-read library method enables short-read second-generation sequencing systems to routinely generate highly accurate and economical long-range sequencing information.* Genome Research 2020. [doi:10.1101/gr.260380.119](https://doi.org/10.1101/gr.260380.119)

Additional source read: *TELL-Seq Library Sequencing User Guide*, Rev B, Universal Sequencing Technology, for custom-primer and run-cycle requirements.

## Evidence record

- 🟢 0.1–5 ng HMW genomic DNA, 3–10 million clonally barcoded magnetic TELL beads and MuA transpososomes react in one unpartitioned tube.
- 🟢 MuA attack forms stable strand-transfer complexes (STCs). Hybridization simultaneously captures those complexes on bead oligos; viscosity and non-barcoded spacer beads restrict diffusion.
- 🟢 A second transpososome tags DNA between STCs to introduce the other library-amplification priming site. STCs are then broken, beads washed, and libraries PCR-amplified off beads while adding P5/P7.
- 🟢 Each 3-µm bead carries about 50,000 barcode templates. The molecular barcode is 18 bases at Index 1 and comes from a >2.4-billion-member pool; the sample index is 8 bases at Index 2.
- 🟢 The defining runs used paired 146-base genomic reads, I1 18 cycles, and I2 8 cycles. Human libraries received 8 PCR cycles; microbial libraries received 13–14.
- 🟢 The vendor guide states that custom R1, R2 and I1 primers are required. I2 needs a custom primer only on some instruments; MiSeq, HiSeq 2000/2500 and NovaSeq use grafted P5 for that read.
- 🔴 Public sources name the proprietary adaptor and custom-primer roles but do not print their sequences. The schematic therefore does not imply standard TruSeq or Nextera primer binding.

## Final construct

The public Figure 1 establishes a paired-end Illumina library with an 18-base TELL molecular barcode at I1 and an 8-base sample index at I2. Reads sharing an I1 barcode are linked back to a long input molecule; the short genomic fragments themselves remain ordinary paired-end inserts.
