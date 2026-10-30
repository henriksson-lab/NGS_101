# TERA-seq — capped source-qualified notes

Ibrahim et al. (2021), TERA-Seq: true end-to-end sequencing of native RNA molecules for transcriptome characterization. DOI10.1093/nar/gkab713.

🟢 Published main Methods59–73/Results147,163 and the complete published SupplementaryTable1 were read. The supplementary PDF identity/Note pp1–2 was inspected as the third original; other figures/Table2 were not reviewed. Mandatory whole-bundle primer scraping completed. XLSX rich-text underlining was checked directly in the original workbook, rather than inferred from flattened extraction.

This is a native-RNA ONT preparation with optional terminal-chemistry selection. Six routes remain branches of one method:

- 5P: poly(A)-selected RNA, on-bead T4RNL1 addition of the5′RNA adapter, washes/water release and RNAClean recovery. This route does **not** subsequently enrich by streptavidin5′adapter capture; both adapted and non-adapted molecules are sequenced.
- 5OH: the corresponding3′phosphorylated adapter and RtcB are used on poly(A) RNA; no decapping step is invented.
- Cap: QuickCIP/washes →RppH decapping →EDTA/washes →5′adapter ligation. Noncapped accessible5P ends are first dephosphorylated.
- Cap&5P: RppH, without the CIP pretreatment, precedes5′adapter ligation.
- TERA3: rRNA/smallRNA subtraction by antisense biotinylated probes; retain the physical unbound fraction →3′adapter ligation →RNAClean →customRTA.
- 5TERA3: rRNA/smallRNA subtraction →RppH/EDTA/RNAClean →5′adapter ligation →MyOneC1 capture/washes →on-bead3′adapter ligation/washes →95%formamide release65°C5min/RNAClean →customRTA. The graph default selects this both-end route.

🟢 Source oligos, 5′→3′:

- 5′adapter: `/5PCBio/rArArUrGrArUrArCrGrGrCrGrArCrCrArCrCrGrArGrArUrCrUrArCrArCrUrCrUrUrUrCrCrCrUrArCrArCrGrArCrGrCrUrCrUrUrCrCrGrArUrCrU`; the5OH variant adds `/3Phos/`. PCBIO is photocleavable biotin. The identical Illumina-like base order is used here as **RNA**; it is not an Illumina sequencing primer or a sample identifier. Code imports the canonical base order, while retaining the locally established all-RNA substrate and end modifications.
- 3′adapter: `5Phos/rGrUrGrUrCrArGrUrCrArCrUrUrCrCrArGrCGG`. Its RNA block precedes two DNA G residues. No ordinary-DNA flattening is allowed by the local constructor.
- RTA-top: `/5Phos/GGCTTCTTCTTGCTCTTAGGTAGTAGGTTC`.
- RTA-custom bottom: `GAGGCGAGCGGTCAATTTTCCTAAGAGCAAGAAGAAGCCCCGCTGGAAGTGACTGACAC`. The final `CCGCTGGAAGTGACTGACAC` is underlined in the original and complements the3′adapter. The local Scene computes the source-known DNA duplex core alignment and retains the unequal free arms; it does not claim that the full RTA strands are mutual reverse complements or reconstruct the motor adapter.

🟢 All routes proceed through RTA ligation by T4DNA ligase →SSIII RT →RNAClean →sequencing-adapter ligation by T4DNA ligase →RNAClean/21µl. RT creates a companion DNA strand: native RNA remains the sequencing substrate. One microliter is consumed for the source Qubit DNA-HS assay; the retained fraction is acquired with MinION/R9.4FLO-MIN106 for48or72h. No PCR, cap-selective perfect purity, photocleavage step or computationally selected physical fraction is invented.

🟡 Historical SQK-RNA001/RNA002, oligo-dT kit and cited antisense-probe components remain named, source-scoped interfaces. Complete enzyme/end chemistry, motor architecture and exact molecule-specific read success remain unassembled. Source short16nt5′adapter is explicitly reference-only/not recommended in TableS1 and is not the model default. No public promotion, accepted absence or inherited whole ONT recipe.

Deferred source scope:

- HistoricalSQK-RNA001/SQK-RNA002 originals and corresponding motor/RTA details; exact product/version applicability unreviewed.
- Published reference9/43 antisense-probe subtraction component originals; source main establishes retained unbound fraction, not full probe chemistry.
- SupplementaryTable2 and supplementaryFigures1–6 remain unreviewed; no absence or independent irrelevance conclusion.
