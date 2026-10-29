# MULTI-seq: sample multiplexing for single-cell RNA sequencing using lipid-tagged indices

Defining DOI: https://doi.org/10.1038/s41592-019-0433-8

## Sources and version boundary

🟢 McGinnis et al., Nature Methods16:619–626(2019), DOI10.1038/s41592-019-0433-8. The primary local workflow follows the PUBLISHED supporting-information protocol in publisherMOESM1, PDFpp33–35. Original renderedp33 was checked against its text layer. The associated PMC6837808 article is explicitly an accepted manuscript, not a verified final-VOR copy; Methods149–158,live-cell labeling169 and library preparation185–229 were read for comparison. The mandatory scraper covered the fetched bundle; available large data workbooks were not chemistry-reviewed.

The sources give distinct oligo versions:

| Component | Published SI protocol,p33 | Associated accepted manuscript |
| --- | --- | --- |
| Anchor5′→3′ | `{Lipid}-TGGAATTCTCGGGTGCCAAGGgtaacgatccagctgtcact` | `{LA/Chol-TEG}-GTAACGATCCAGCTGTCACTTGGAATTCTCGGGTGCCAAGG` |
| Additive primer5′→3′ | `CTTGGCACCCGAGAATTCC` | `CCTTGGCACCCGAGAATTCC` |

These differences reverse anchor domain order and add/remove one5′C from the additive primer. They are SOURCE-VERSION differences, not an established conflict between two final published originals or a publication-error finding. The graph consistently uses the published SI sequence/procedure. A separately named accepted-manuscript interface is shown only for comparison; final-VOR equivalence and which variant was used in a specific dataset remain unresolved. No recipe is merged across versions. SI uses generic lipid labels; the selected model does not fill their exact chemical identity from the manuscript's LA/PA description.

## Physical oligo complex and identifiers

🟢 Co-anchor: `AGTGACAGCTGGATCGTTAC` with a source3′lipid. Barcode: `CCTTGGCACCCGAGAATTCCANNNNNNNNA30`. Scene computes pairing of separate anchor,co-anchor and barcode strands for each source version; these do not become one covalent DNA strand. Barcode polyA and sample-tag regions are unpaired overhangs in the partial view, not extra complementarity to an overlapping domain.

The8base membrane barcode identifies the biological sample. It is separate from the later6base TruSeqRPIX sequencer library index. Cell barcode and UMI derive from the source capture workflow and remain unassigned-width in this partial model; no arbitrary current10x architecture is inserted.

## Published SI route

🟢 Mix anchor/barcode1:1 before cell labeling. Associate5min on ice without BSA/serum, add co-anchor and wait5min, then add ice-cold1%BSA/PBS to quench free lipidoligos. Quenching is not physical disposal: BSA/reagent species persist until source centrifugation/washes. Wash at least twice, combine/filter samples and follow the historical10x workflow through cDNA retrieval. This handoff is explicitly unassembled. Synthetic DNA barcode molecules are not renamed RNA merely because the capture process also handles transcripts.

🟢 Amplification adds the SI19nt primer to the historical10xcDNA primer mix. At0.6×SPRI separation, endogenouscDNA remains on beads and barcode-derivedDNA remains in the SUPERNATANT; both fractions are retained. The endogenous branch remains a separate unresolved10x library handoff. Recover barcodeDNA using3.2×finalSPRI and1.8×isopropanol, then ethanol washes/EBelution.

🟢 Use3.5ng barcodeDNA in KAPAHiFiPCR with SIUniversalI5/P5 and TruSeqRPIX. The SI prints an instruction to repeat98/60/72°C steps8–12times; the model preserves that instruction rather than asserting a total-cycle interpretation. RPIX is `CAAGCAGAAGACGGCATACGAGATNNNNNNGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA`. Final1.6×SPRI removes primers/contaminants physically, without assuming complete purity.

## Deferrals and lineage

Full historical10xv2 capture/RT/cDNA architecture, source barcode whitelist/run formats and final read-primer landing remain unassembled. Published main VOR and study-specific reconciliation of source variants remain unresolved; no absence claim is made. Main Methods explicitly bases its capture/library workflow on10xV2 and CITE-seq; SI provides its own workflow. The sample-tag architecture and Weber membrane-anchor strategy are component ancestry suggestions, not authority to reuse an antibody-conjugation or entire10x recipe. CMO and nuclei variants are separate unmodeled routes.

## Completion

Bounded original-source pass only. Catalogue remains notes. Private partial graphs retain explicit uncovered source inventory and unresolved terminal handoffs; passing structural/local checks does not complete source review. No two-reader missing-information pass has been performed.
