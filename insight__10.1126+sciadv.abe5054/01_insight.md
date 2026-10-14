# INSIGHT — research notes

Defining paper: Wu *et al.* (2021), DOI [10.1126/sciadv.abe5054](https://doi.org/10.1126/sciadv.abe5054). The PMC full text and protocol table were read.

## Evidence

- 🟢 Table 2 gives the P8 target primers, 5-nt sample barcodes, T7 promoter, Illumina handle, outer indexed PCR primers and custom T7-containing sequencing primer.
- 🟢 NASBA uses reverse transcriptase, RNase H and T7 RNA polymerase at 41 °C.
- 🟢 Completed reactions are pooled only after the sample barcode pair is installed; residual primers are removed before a minimal-cycle universal PCR.
- 🟢 The same Table 2 prints the alternative readout oligos, which do not enter the sequenced library: a 2′-O-methyl toehold beacon `FAM-AUUGACAGUCUACUAAUUUGGUUAAAAACAAAUGUGUCAA-BHQ1dT-UUCAACUUCAAUG-propyl`, PCRD capture probe A `FAM-AAAAGUCUACUAAUUUGGUUAAAA`, and probe B `ACAAAUGUGUCAAUUUCAACUUCA-Biotin`. These are RNA sequences as written; the beacon's internal quencher and terminal propyl are retained verbatim.

The lateral-flow and molecular-beacon readouts are useful parallel endpoints, but the schematic concentrates on the molecules that enter NGS.
