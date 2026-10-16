# scSNV-seq — source notes

Defining paper: Cooper et al., *Genome Biology* (2024), DOI
[10.1186/s13059-024-03169-y](https://doi.org/10.1186/s13059-024-03169-y).

🟢 The vector carries a 12-nt transcribed PuroR barcode and a 6-nt iBAR in the guide.
The workflow combines targeted single-cell DNA genotyping with 10x 5′ transcriptomes;
guide/iBAR and PuroR libraries use spiked-in RT primers and separate nested enrichments.

## Oligo audit

🟢 Sources checked: the defining PMC article and JATS, associated bioRxiv manuscript,
all five supplementary files, and especially Supplementary Table S2. That table publishes
the vector-construction oligos, the complete JAK1 Tapestri amplicon panel, both custom 10x
RT spike-ins, and the shared indexed nested primer.

🟢 The guide spike-in is oJR160,
`AAGCAGTGGTATCAACGCAGAGTACCAAGTTGATAACGGACTAGCC`; the PuroR spike-in is
`AAGCAGTGGTATCAACGCAGAGTACggcttgtactcggtcatggt`. Their 5′ SMART handle is
followed by a scaffold- or PuroR-specific binding region.

🟢 `Kosuke_Guide_i5longFwd` and `Puro_i5longFwd` are the same published primer:
`AATGATACGGCGACCACCGAGATCTACACNNNNNNNNNNACACTCTTTCCCTACACGACGCTCTTCCGATCT`,
where the ten `N` bases are the library index.

🟡 The paper does not print the paired commercial 10x 5′ primer/bead construct or the
complete Mission Bio Tapestri sequencing shell. The exact custom oligos and target-panel
primers are shown separately. Historical vendor guides and the authors' linked code
still need a second source check before those commercial endpoints can be called absent.
