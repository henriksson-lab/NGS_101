# DBiT-seq

## 1. What it is

DBiT-seq installs two spatial barcode coordinates directly in a fixed tissue section using perpendicular microfluidic flows, then reads RNA and optionally antibody-derived protein tags by NGS.

## 2. Sources and evidence

- 🟢 Liu Y et al. *High-Spatial-Resolution Multi-Omics Sequencing via Deterministic Barcoding in Tissue.* Cell 183:1665–1681 (2020). [doi:10.1016/j.cell.2020.10.026](https://doi.org/10.1016/j.cell.2020.10.026).
- 🟢 Liu Y et al. *Deterministic barcoding in tissue for spatial omics sequencing (DBiT-seq).* STAR Protocols 2:100532 (2021). [doi:10.1016/j.xpro.2021.100532](https://doi.org/10.1016/j.xpro.2021.100532). Supplementary Table S2 prints the complete oligos used below.
- 🟡 The defining 2020 article links its own Table S2 oligo workbook. The local PMC author manuscript contains the article text but that workbook is behind PMC's supplement download gate. The 2021 author table is direct evidence for the sequences below; exact identity of the 2020 and 2021 tables is still under review. The linked 2019 bioRxiv preprint also says its exact oligos are in a separate Excel Table S2, not in its main PDF.
- 🟢 Flow A uses 50 oligo-dT reagents, each with an 8-nt A barcode and 15-nt linker, during in situ reverse transcription.
- 🟢 Perpendicular flow B supplies a 15-nt linker, 8-nt B barcode, 10-nt UMI, 22-nt PCR handle and 5′ biotin. T4 ligation creates A/B pixel addresses.
- 🟢 Barcode A-1 is `/5Phos/AGGCCAGAGCATTCGAACGTGATTTTTTTTTTTTTTTTVN`; Barcode B-1 is `/5Biosg/CAAGCGTTGGCTTCTCGCATCTNNNNNNNNNNAACGTGATATCCACGTGCTTGAG`.
- 🟢 The universal ligation linker is `CGAATGCTCTGGCCTCTCAAGCACGTGGAT`. Its 5′ half is the reverse complement of the A 5′ linker (`AGGCCAGAGCATTCG`), and its 3′ half is the reverse complement of the B 3′ linker (`ATCCACGTGCTTGAG`). Equivalently, the reverse complement of the complete linker is B-linker followed by A-linker. It therefore anneals the B 3′-OH beside the A 5′ phosphate and exposes one central nick for T4 DNA ligase.
- 🟢 Recovered cDNA is template-switched, PCR-amplified, Nextera-XT-tagmented and sequenced 2×100.
- 🟢 The same STAR Protocols Table S2 "PCR Oligos" sheet also prints the template-switch oligo `AAGCAGTGGTATCAACGCAGAGTGAATrGrG+G` (terminal locked G), Primer 1 `CAAGCGTTGGCTTCTCGCATCT`, and Primer 2 `AAGCAGTGGTATCAACGCAGAGT`. Its optional antibody-tag branch lists Primer 3 `CCTTGGCACCCGAGAATT*C*C` (two phosphorothioate linkages). These modifications should not be lost when ordering the assay reagents.
- 🟢 The table includes complete P5 and indexed P7 oligos for the author's detailed workflow; for example P5 `AATGATACGGCGACCACCGAGATCTACACTAGATCGCTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG` and P7-701 `CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGCAAGCGTTGGCTTCTCGCATCT`.
- 🟡 Barcode A-1 and B-1 demonstrate the geometry. The other channel-specific whitelist members remain in Supplementary Table S2 rather than being recopied here.

## 3. Molecular path

The first barcode is covalently installed by RT in one channel direction. The second is installed by ligation in the perpendicular direction. Consequently the A/B combination is a physical pixel address rather than an inferred coordinate.

The two barcode oligos do not ligate through unspecified generic ends. Their published
15-nt terminal arms form a continuous 30-bp duplex with the universal splint. Their long
barcode/UMI/primer bodies project beyond opposite ends of that splint; only the central
B–A nick is sealed.
