# DBiT-seq

## 1. What it is

DBiT-seq installs two spatial barcode coordinates directly in a fixed tissue section using perpendicular microfluidic flows, then reads RNA and optionally antibody-derived protein tags by NGS.

## 2. Sources and evidence

- 🟢 Liu Y et al. *High-Spatial-Resolution Multi-Omics Sequencing via Deterministic Barcoding in Tissue.* Cell 183:1665–1681 (2020). [doi:10.1016/j.cell.2020.10.026](https://doi.org/10.1016/j.cell.2020.10.026).
- 🟢 Su G et al. *Spatial multi-omics sequencing for fixed tissue via DBiT-seq.* STAR Protocols 2:100532 (2021). [doi:10.1016/j.xpro.2021.100532](https://doi.org/10.1016/j.xpro.2021.100532). Supplementary Table S2 prints the complete oligos used below.
- 🟡 The defining 2020 article links its own Table S2 oligo workbook. The local PMC author manuscript contains the article text but that workbook is behind PMC's supplement download gate. The 2021 author table is direct evidence for the sequences below; exact identity of the 2020 and 2021 tables is still under review. The linked 2019 bioRxiv preprint also says its exact oligos are in a separate Excel Table S2, not in its main PDF.
- 🟢 Flow A uses 50 oligo-dT reagents, each with an 8-nt A barcode and 15-nt linker, during in situ reverse transcription.
- 🟢 Perpendicular flow B supplies a 15-nt linker, 8-nt B barcode, 10-nt UMI, 22-nt PCR handle and 5′ biotin. T4 ligation creates A/B pixel addresses.
- 🟢 Barcode A-1 is `/5Phos/AGGCCAGAGCATTCGAACGTGATTTTTTTTTTTTTTTTVN`; Barcode B-1 is `/5Biosg/CAAGCGTTGGCTTCTCGCATCTNNNNNNNNNNAACGTGATATCCACGTGCTTGAG`.
- 🟢 The universal ligation linker is `CGAATGCTCTGGCCTCTCAAGCACGTGGAT`. Its 5′ half is the reverse complement of the A 5′ linker (`AGGCCAGAGCATTCG`), and its 3′ half is the reverse complement of the B 3′ linker (`ATCCACGTGCTTGAG`). Equivalently, the reverse complement of the complete linker is B-linker followed by A-linker. It therefore anneals the B 3′-OH beside the A 5′ phosphate and exposes one central nick for T4 DNA ligase.
- 🟢 Recovered cDNA is template-switched, PCR-amplified, Nextera-XT-tagmented and sequenced 2×100.
- 🟢 The same STAR Protocols Table S2 "PCR Oligos" sheet also prints the template-switch oligo `AAGCAGTGGTATCAACGCAGAGTGAATrGrG+G` (terminal locked G), Primer 1 `CAAGCGTTGGCTTCTCGCATCT`, and Primer 2 `AAGCAGTGGTATCAACGCAGAGT`. Its optional antibody-tag branch lists Primer 3 `CCTTGGCACCCGAGAATT*C*C` (two phosphorothioate linkages). These modifications should not be lost when ordering the assay reagents.
- 🟢 The table includes complete P5 and indexed P7 oligos for the author's detailed workflow; for example P5 `AATGATACGGCGACCACCGAGATCTACACTAGATCGCTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG` and P7-701 `CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGCAAGCGTTGGCTTCTCGCATCT`.
- 🟢 STAR Protocols steps 60–66: 0.75–1 ng cDNA, Tagment DNA buffer + Amplicon Tagment Mix 55 °C 5 min, NT buffer, then 12-cycle PCR with the P5 primer and an indexed P7 primer, 0.7× KAPA Pure cleanup. Every indexed P7 oligo (701–706) ends `…AGATGTGTATAAGAGACAG` + `CAAGCGTTGGCTTCTCGCATCT`, i.e. its 3′ end is the Barcode B handle; the P5 oligo ends in the Nextera S5/mosaic end and carries a fixed i5 `TAGATCGC`.
- 🟢 Cell STAR Methods "Template switch and PCR amplification": Maxima H Minus TSO reaction on streptavidin C1 beads (RT 30 min, 42 °C 90 min), Kapa HiFi HotStart 5 cycles on beads, then qPCR-monitored PCR to plateau; 0.6× AMPure keeps mRNA cDNA (>300 bp). "Sequencing library was then built with a Nextera XT kit … HiSeq 4000 … pair-end 100×100". "Quantification": Read 2 is processed to extract UMI, Barcode A and Barcode B; Read 1 is mapped.
- 🟡 Derived orientation (computed by the model): the B-handle end is the P7/Read 2 end, so Read 2 crosses the 22-nt handle, then the UMI, Barcode B, the two linkers and Barcode A, consistent with the Cell analysis. The first-strand 3′ end carries the reverse complement of the template-switch oligo, which Primer 2 primes.
- Model correction (2026-10-17 bounded pass): the page previously drew the spatial address at the Read 1 end and appended the Primer 2 sequence (not its complement) to the first strand. Both were repository modelling errors, not source discrepancies.
- Not drawn: the ADT branch (0.6× supernatant → 2× AMPure → 10-cycle PCR with P7 and "P5 cite" primers, pooled 1:9 with the mRNA library; Cell STAR Methods), Proteinase K lysis detail is a graph condition only.
- 🟡 Barcode A-1 and B-1 demonstrate the geometry. The other channel-specific whitelist members remain in Supplementary Table S2 rather than being recopied here.

## 3. Sources inspected in the 2026-10-17 bounded pass

- Cell 2020 author manuscript PMC7736559 JATS: Results line 98 and 106, STAR Methods lines 194, 214–236 (ADT staining, flow A RT, flow B ligation, lysis/streptavidin, TSO/PCR, amplicon purification/Nextera XT/HiSeq), 252–254 (read processing).
- STAR Protocols 2021 PMC8132129 JATS steps 60–67 (library preparation) and Table S2 workbook (mmc2: PCR Oligos, Barcode A, Barcode B sheets).
- Deferred: Cell Table S2/S3 workbooks (NIHMS1645400-supplement-10.xlsx and others, PMC download gate; manual URLs in the source manifest), so 2020-vs-2021 oligo identity and the 2020 index primers remain open.

## 4. Molecular path

The first barcode is covalently installed by RT in one channel direction. The second is installed by ligation in the perpendicular direction. Consequently the A/B combination is a physical pixel address rather than an inferred coordinate.

The two barcode oligos do not ligate through unspecified generic ends. Their published
15-nt terminal arms form a continuous 30-bp duplex with the universal splint. Their long
barcode/UMI/primer bodies project beyond opposite ends of that splint; only the central
B–A nick is sealed.

## Native ownership migration, 2026-10-18

This pass read three originals: the defining Cell 2020 article (PMC7736559 JATS Results 98, STAR Methods 214–236 and read processing 252–254), the later STAR Protocols 2021 article (PMC8132129 detailed cDNA cleanup 45–59 and library preparation 60–67), and that later article's original Table S2 workbook (PCR Oligos and A/B barcode sheets). The source manifest assigns the cached mmc2 workbook to the 2021 implementation; its publisher-prefixed duplicate filename does not establish identity with Cell 2020 Table S2. The original 2020 workbook remains a deferred source target. Ran the full bundle oligo scraper and verified the used oligo rows in the actual 2021 workbook.

The detailed graph is explicitly the **2021 implementation reference**, retaining the 2020 paper as the defining method; its historical PE100/100 run is documented separately rather than transferred to the later implementation. It does not establish 2020/2021 oligo-table identity. In particular, the graph uses the 2021 0.8× cDNA cleanup; Cell 2020 uses a 0.6× modality fractionation. These are version-specific recipes, not a publication discrepancy or interchangeable ratios. The optional 2020 antibody-tag branch remains outside this mRNA reference graph.

The local endpoint comes directly from the 2021 P5/P7-701 oligos and A-1/B-1 chain: fixed P5-side index TAGATCGC and representative P7-701 index are explicit features. No generic Nextera endpoint selects index lengths, values or fragment ends. A new complete A-1 transcription check found that its literal poly(dT) tract has **15 T bases after the eight-base barcode**, whereas the repository had used 16. The 2021 reference is corrected; this is a repository transcription error and does not assign the unrecovered 2020 oligo.

RT output and the template-switched product are single strands, distinct from the amplified duplex. Both PCR geometries now use their actual incoming precursor Scenes. Native first-cycle/tagmented-end geometry remains unassigned; a finished library is never substituted as the template. The graph adds the source cDNA and final 0.7× KAPA Pure cleanups. Candidate fragments lacking the B handle are unresolved material, since failure to amplify does not establish physical removal. Acquisition is unspecified for this 2021 geometry reference, including both reads and all cycle counts; no 2020 run values generate cycle spans on its endpoint. Curator-written reviewed targets remain independent of graph enumeration. No absence conclusion is accepted.

### Independent version/run correction

Rechecking the existing Su et al. 2021 original's header and library steps 60–67 identified the correct first author and clarified that the inspected later implementation passage establishes library construction/QC, without positively assigning a run recipe. The 2021 graph and rendered sequencing section therefore use acquisition-unspecified geometry with no read-length values. Liu et al. Cell 2020's HiSeq4000 PE100/100 remains a separate historical statement. This correction adds no new original beyond the three already inspected for the migration and does not infer absence from those passages.
