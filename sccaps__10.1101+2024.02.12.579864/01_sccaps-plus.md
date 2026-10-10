# scCAPS+

## Sources read

- 🟢 Chen et al., bioRxiv preprint [doi:10.1101/2024.02.12.579864](https://doi.org/10.1101/2024.02.12.579864), the shared defining record.
- 🟢 Chen et al., version of record, *Genome Biology* 26:244 (2025), [doi:10.1186/s13059-025-03708-1](https://doi.org/10.1186/s13059-025-03708-1), especially Methods and Additional file 5.

## Molecular path

🟢 Single cells are FACS-sorted, lysed and tagmented with one of 96 preassembled i5/i7-barcoded Tn5 complexes. After SDS stripping, neutralization and 72 °C gap filling, 96 cell-barcoded reactions are pooled.

🟢 ACT+BF4− selectively oxidizes 5hmC to 5fC. Pinnick oxidation with sodium chlorite and 2-methyl-2-butene converts 5fC to 5caC. Pyridine–borane reduces 5caC to DHU, which PCR reads as T. 5mC and unmodified C remain C, making scCAPS+ selective for 5hmC.

🟢 Nextera XT index primers amplify the converted pool. The paper uses 120-bp paired-end sequencing with custom NovaSeq primers; Additional file 5 Tables S3–S4 define the transposome and sequencing oligos.

## Protocol boundary

This is a separate schematic from scTAPS despite the joint paper. The cell-barcoding front end is shared, but scCAPS+ has two selective chemical oxidations and a 5hmC-only readout; scTAPS instead uses mTET1 and reports 5mC+5hmC.
