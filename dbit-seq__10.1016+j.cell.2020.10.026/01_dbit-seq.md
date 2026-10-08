# DBiT-seq

## 1. What it is

DBiT-seq installs two spatial barcode coordinates directly in a fixed tissue section using perpendicular microfluidic flows, then reads RNA and optionally antibody-derived protein tags by NGS.

## 2. Sources and evidence

- 🟢 Liu Y et al. *High-Spatial-Resolution Multi-Omics Sequencing via Deterministic Barcoding in Tissue.* Cell 183:1665–1681 (2020). [doi:10.1016/j.cell.2020.10.026](https://doi.org/10.1016/j.cell.2020.10.026).
- 🟢 Flow A uses 50 oligo-dT reagents, each with an 8-nt A barcode and 15-nt linker, during in situ reverse transcription.
- 🟢 Perpendicular flow B supplies a 15-nt linker, 8-nt B barcode, UMI, 22-nt PCR handle and biotin. T4 ligation creates A/B pixel addresses.
- 🟢 Recovered cDNA is template-switched, PCR-amplified, Nextera-XT-tagmented and sequenced 2×100.
- 🟡 The UMI is drawn as a role placeholder; whitelist sequences remain in the paper’s supplementary table rather than being recopied here.

## 3. Molecular path

The first barcode is covalently installed by RT in one channel direction. The second is installed by ligation in the perpendicular direction. Consequently the A/B combination is a physical pixel address rather than an inferred coordinate.
