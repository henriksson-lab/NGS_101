# scTAPS

## Sources read

- 🟢 Chen et al., bioRxiv preprint [doi:10.1101/2024.02.12.579864](https://doi.org/10.1101/2024.02.12.579864), the catalogue's defining record.
- 🟢 Chen et al., version of record, *Genome Biology* 26:244 (2025), [doi:10.1186/s13059-025-03708-1](https://doi.org/10.1186/s13059-025-03708-1), especially Methods and Additional file 5.

## Molecular path

🟢 Single cells are FACS-sorted into 96 wells, lysed and tagmented at 50 °C for 15 min with one of 96 preassembled barcoded Tn5 complexes. SDS strips Tn5; NP-40 neutralizes SDS; a 72 °C extension fills the tagmentation gaps. The 96 reactions are then pooled.

🟢 Two mTET1 oxidation rounds convert 5mC and 5hmC to 5caC. Pyridine–borane reduces 5caC to DHU; uracil-tolerant PCR copies DHU as T while unmodified C remains C. Thus scTAPS directly reports the combined 5mC+5hmC fraction.

🟢 Nextera XT index primers amplify the converted pool. The paper uses 120-bp paired-end sequencing and custom NovaSeq sequencing primers; Additional file 5 Tables S3–S4 are authoritative for the 96 transposome oligos and custom primers.

## Oligo source review

🟢 Additional file 5, Table S3 prints eight i5 and twelve i7 Tn5 transfer oligos,
which form 96 pairwise cell labels. Every listed transfer oligo begins `/5Phos/`;
the shared reverse mosaic-end strand is `/5Phos/CTGTCTCTTATACACATCT`.
The first i5 member is
`/5Phos/TCGTCGGCAGCGTCTCCACGCTATAGCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG`;
the first i7 member is
`/5Phos/GTCTCGTGGGCTCGGCTGTCCCTGTCCGAGTAATCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG`.
These are ordered transposome oligos, not the Nextera XT PCR primers.

🟢 Table S4 gives four custom sequencing primers: Read 1
`GCGATCGAGGACGGCAGATGTGTATAAGAGACAG`, Read 2
`CACCGTCTCCGCCTCAGATGTGTATAAGAGACAG`, Index 1
`CTGTCTCTTATACACATCTGAGGCGGAGACGGTG`, and Index 2
`CTGTCTCTTATACACATCTGCCGTCCTCGATCGC`.
The three oligos in Table S1 make a 5hmC spike-in control; they are not genomic
library adapters. 🟡 The current schematic uses a generic Nextera completion
shell and therefore does not yet reconstruct the published custom primer sites.

## Protocol boundary

scTAPS and scCAPS+ share one paper and front-end barcoding, but they are separate protocols: scTAPS uses enzymatic TET oxidation and reports 5mC+5hmC, whereas scCAPS+ uses selective chemical oxidations and reports 5hmC alone.
