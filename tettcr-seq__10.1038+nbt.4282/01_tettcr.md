# TetTCR-seq

Defining original: **High-throughput determination of the antigen specificities of T cell receptors in single cells**, DOI[10.1038/nbt.4282](https://doi.org/10.1038/nbt.4282).

Capped source/model pass,2026-10-29. Private and incomplete;no accepted source-absence claim. Structural graph checks are separate from chemistry/source completeness.

## Selected originals

- `Zhang_Nat_Biotech_2018.pdf`: Institutional mirror printed title,authors,doi10.1038/nbt.4282 and12Nov2018 publication/December2018journal identity verified before chemistry. Selected OnlineMethods973–1064: linker/protein,dual DNA-template uses,UVpMHC exchange,tetramer/cell handling,retained flowthrough and jointRT/PCR then parallel TCR/DNA libraries.
- `nbt.4282_supp_NIHMS1508189-supplement-11.xlsx`: Published MOESM15 Table10 has a manifest-validated byte-identical alias to this cached file. Selected source DNA Linker manufacturing syntax/5AmMC12/iSp18/12N,HCV-KLV(WT) template,P1/P2,V1F/V1Rmember1 andILLU_f/r. This alias establishes selected file identity,not equivalence of all author-manuscript sources.
- `nbt.4282_41587_2018_BFnbt4282_MOESM4_ESM.pdf`: Published SI Fig3 caption125–130 confirms covalent linker/fluorescentstreptavidin and DNA-template extension;control classification48–65 distinguishes tetramernegative single cells from enriched fraction and flowthrough fluorescence background. Remaining figures/tabledata not fully audited.

## Chemistry and limits

🟢 The university mirror is a published NatureBiotechnology original with printed title,DOI and date verified;the accepted PMC manuscript is not used to replace it. Table10's publisherMOESM15 download is byte-identical to the cached NIHMS supplement11,as recorded in the canonical manifest. That is exact file provenance,not a claim of whole-version equivalence.

🟢 The peptide-encoding DNA has two physically separate uses. One aliquot is PCR-amplified for PURExpress IVTT,which retains a mixture of DNA template,RNA,peptide,enzymes and reagents. A second aliquot encodes its corresponding fluorescentstreptavidin DNA linker by annealing55°C then exo-minusKlenow37°C1h and EDTA quenching. No purification is invented at EDTA addition. The local graph preserves these independent preparation branches before matching the source antigen and DNA reporter.

🟢 Source DNA Linker as written is `/5AmMC12/ /iSp18/ TAG TAC TCA GAG GTT GAT CTA CAT TG (N:25252525)(N)(N) (N)(N)(N) (N)(N)(N) (N)(N)(N) GAC GAT GAC GAC AAG`. There are12explicitNsymbols;the first balanced-synthesis directive does not denote25252525bases. The5′amineC12 andiSp18 are nonnucleotide attachment/spacer chemistry,not12or18 extra DNA bases. The source amine is not asserted free after conjugation. Synthetic linkerMID12 counts reporter molecules;it precedes PCR and is distinct from well identity,TCR RNA molecule identity and laterNHNHN windows.

🟢 PE/APCstreptavidin–linker conjugates are cleaned by six100kDaVivaspin washes,concentration and0.2µm filtration. Solulink attachment geometry remains unread. Source biotinylated UV-labilepMHC is added directly to the completed IVTT mixture,UVexchanged60min on ice then incubated4°C≥12h. This is peptide exchange,not DNA photorelease. Matched DNA-streptavidin and pMHC are incubated1:6.7,4°C1h;that input ratio is not a predicted exact output tetramer stoichiometry. Panel pooling/100kDa concentration is a physical operation;remaining species and recovery are unquantified.

🟢 SourceCD8 cells are stained~60nM each tetramer4°C1h,washed,incubated with antiPE/antiAPCmicrobeads and selected byLScolumn. The enriched and tetramer-depleted flowthrough fractions are both retained. Flowthrough phenotyping is not assigned an NGS library. Published SI uses tetramernegative single cells from the enriched fraction for MID-background classification;binding positivity is not predicted from reporter presence. One representative sorted cell enters4µl lysis,while other cell/sort fractions retain unresolved fates.

🟢 TCRRT occurs in the same lysis well with DNA-BCs present. FirstPCR includesP1/P2at100nM for joint TCR/DNA-BC amplification. Subsequent preparations use separate parallel wells,not physically joined TCR–antigen amplicons. DNA-BCsecondPCR uses1µl firstPCR in5µl total,ExTaqHS,10cycles95/55/72 with V1F/V1R well tags. TCR aliquot and unallocated firstPCR remainder remain explicit. Third DNA-BCPCR usesILLU_f/ILLU_r for5cycles;multiple wells are pooled and gel purified. TCR continues through the cited parent PCR method,then pooling/gel and MiSeqV2. The source DNA-BC depth is≥6000reads/cell;exact read/index cycles and whole-library geometry remain unassembled.

🟢 Table10 source P1 `TCAGAGGTTGATCTACATTG`,P2 `AGCGAGGTGCTTCGTTA`. The selected HCV-KLV(WT)DNA is `CGAGGTGCTTCGTTAAACAGCGTTGATACCCAGAGCAACCAGTTTCTTGTCGTCATCGTC`. V1F/V1Rmember1 inline well tags are `ATCACG` and `TCAAG`;their experimental well assignment is not inferred. The printed5-baseNHNHN patterns retain constrained synthesis notation with unassigned function;they are not called another UMI. ILLU_r prints fixedAA followed by sixX;no full acquired index boundary or8cycle read is inferred solely from that local window.

Component ancestry: the published main explicitly follows prior single-cell TCR protocol(ref7) with local P1/P2 additions and parallel DNA-BC handling. Exact parent recipe/version remains unread before code inheritance. pMHC UVexchange cites earlier protocols5/6;that does not make the source IVTT/tag construction identical to older bulk dextramer library chemistry. Antigen binding,RNA RT,DNA template encoding and later computational association remain distinct.

## Deferred

- Cited single-cell TCR protocol ref7 and complete historical RT/nestedPCR/index/run geometry
- Solulink conjugation guide and pMHC UVexchange predecessor sources5/6
- Complete source peptide/linker/copy duplex geometry,manufacturing annotations and NHNHN purpose
- Remaining published SI/linked originals and independent source review

## Validation

Local source transcription and physical branch/state tests are recorded in the source-review JSON. No finished library or full source review is claimed.
