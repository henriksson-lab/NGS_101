# DBiT-seq

## 1. What it is

DBiT-seq installs two spatial barcode coordinates directly in a fixed tissue section using perpendicular microfluidic flows, then reads RNA and optionally antibody-derived protein tags by NGS.

## 2. Sources and evidence

- 🟢 Liu Y et al. *High-Spatial-Resolution Multi-Omics Sequencing via Deterministic Barcoding in Tissue.* Cell 183:1665–1681 (2020). [doi:10.1016/j.cell.2020.10.026](https://doi.org/10.1016/j.cell.2020.10.026).
- 🟢 Liu Y et al. *Deterministic barcoding in tissue for spatial omics sequencing (DBiT-seq).* STAR Protocols 2:100532 (2021). [doi:10.1016/j.xpro.2021.100532](https://doi.org/10.1016/j.xpro.2021.100532). Supplementary Table S2 prints the complete oligos used below.
- 🟢 Flow A uses 50 oligo-dT reagents, each with an 8-nt A barcode and 15-nt linker, during in situ reverse transcription.
- 🟢 Perpendicular flow B supplies a 15-nt linker, 8-nt B barcode, 10-nt UMI, 22-nt PCR handle and 5′ biotin. T4 ligation creates A/B pixel addresses.
- 🟢 Barcode A-1 is `/5Phos/AGGCCAGAGCATTCGAACGTGATTTTTTTTTTTTTTTTVN`; Barcode B-1 is `/5Biosg/CAAGCGTTGGCTTCTCGCATCTNNNNNNNNNNAACGTGATATCCACGTGCTTGAG`.
- 🟢 The universal ligation linker is `CGAATGCTCTGGCCTCTCAAGCACGTGGAT`. It is the reverse complement of the B 3′ linker (`ATCCACGTGCTTGAG`) followed by the A 5′ linker (`AGGCCAGAGCATTCG`). Thus it anneals the B 3′-OH beside the A 5′ phosphate and exposes one central nick for T4 DNA ligase.
- 🟢 Recovered cDNA is template-switched, PCR-amplified, Nextera-XT-tagmented and sequenced 2×100.
- 🟡 Barcode A-1 and B-1 demonstrate the geometry. The other channel-specific whitelist members remain in Supplementary Table S2 rather than being recopied here.

## 3. Molecular path

The first barcode is covalently installed by RT in one channel direction. The second is installed by ligation in the perpendicular direction. Consequently the A/B combination is a physical pixel address rather than an inferred coordinate.

The two barcode oligos do not ligate through unspecified generic ends. Their published
15-nt terminal arms form a continuous 30-bp duplex with the universal splint. Their long
barcode/UMI/primer bodies project beyond opposite ends of that splint; only the central
B–A nick is sealed.
