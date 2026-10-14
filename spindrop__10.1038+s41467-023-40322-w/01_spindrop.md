# spinDrop

**Defining source:** doi:10.1038/s41467-023-40322-w.

## Molecular boundary

🟢 Cells or nuclei are co-encapsulated with inDrop-v3 microgels, fluorescence-sorted, heat-treated and sent through a second device that picoinjects RT mix.

## Evidence boundary

🟢 Sorting and picoinjection are explicit physical steps; the final identifiers follow the inDrop-v3 molecular scheme.

## Oligo source review

🟢 The linked preprint (doi:10.1101/2023.01.12.523500, extracted methods
around lines 1083–1139) says that the barcoded beads were prepared according
to the cited inDrop protocol with **inDrop v3 oligonucleotide barcoding**.
The sequencing recipe is Read 1 = 61 cycles, Index 1 and 2 = 8 cycles each,
and Read 2 = 14 cycles. The source describes release of bead oligos by UV
exposure and RT after picoinjection, but prints no ordering-level capture
oligo, amplification primer, index primer, or sequencing primer bases.

🟡 The spinDrop molecular model is currently role-level. To resolve its
primer-binding geometry, the cited inDrop preparation protocol (Zilionis et
al., doi:10.1038/nprot.2016.154, ref. 28), the cited v3 barcoding source
(Briggs et al., doi:10.1126/science.aar5780, ref. 29), and the exact bead
implementation used here must be compared; an inDrop sequence from another
version cannot establish the bases of these beads. The article's available
workbook supplements are cell counts, device fabrication and cost tables,
not an oligo table.
