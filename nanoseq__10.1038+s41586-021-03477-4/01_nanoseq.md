# NanoSeq — research notes

Defining source: Abascal *et al.*, [10.1038/s41586-021-03477-4](https://doi.org/10.1038/s41586-021-03477-4).

## Evidence used

- 🟢 NanoSeq is a duplex-sequencing protocol based on restriction-enzyme fragmentation and restrictive end repair.
- 🟢 Dideoxynucleotides prevent extension from pre-existing nicks during end repair, avoiding fixation of damage as a mutation on both descendants.
- 🟢 Duplex adapters associate Watson- and Crick-derived read families before consensus calling.
- 🟢 The Methods specify HpyCH4V fragmentation, A-tailing with dATP plus ddBTPs, and ligation of the named xGen Duplex Seq Adapter (IDT 1080799). The amplified libraries use UDI-containing PCR primers.
- 🟢 Two qPCR oligos are printed in the Methods: `NanoqPCR1` = `ACACTCTTTCCCTACACGAC` and `NanoqPCR2` = `GTGACTGGAGTTCAGACGTG`. They measure adapter-ligated library concentration; they are not the UDI-containing library PCR primers.
- 🟢 The read-processing Methods specify a three-base random barcode at the beginning of each read, followed by the adapter overhang and HpyCH4V restriction-end bases (`NNNTCA`, with an optional spacer base before `TCA`). Those observed bases constrain the library end but do not specify the complete adapter duplex or its chemical modifications.
- 🔴 Neither the paper's Methods nor its Supplementary Notes/Tables print the complete xGen adapter strands or the particular UDI-containing amplification primers. The article points to an author Protocol Exchange method (Lensing *et al.*, [10.21203/rs.3.pex-1298/v1](https://doi.org/10.21203/rs.3.pex-1298/v1)); until its full oligo appendix can be checked, the end-to-end NanoSeq library remains source-limited.

## Protocol boundary

NanoSeq is not Nanopore sequencing. It is a distinct Illumina rare-variant protocol and does not replace the general Duplex Sequencing page.
