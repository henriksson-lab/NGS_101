# Perturb-Seq: Dissecting Molecular Circuits with Scalable Single-Cell RNA Profiling of Pooled Genetic Screens

Defining DOI: https://doi.org/10.1016/j.cell.2016.11.038

## Identity and defining distinction

🟢 Dixit et al., Cell167:1853–1866.e17(2016), DOI10.1016/j.cell.2016.11.038. Source is the associated final accepted manuscript PMC5181115, explicitly marked author-manuscript in the external manifest; final-VOR equivalence is not asserted. Read Results105–113, vector construction256–316, single-cell library preparation318–379, and dictionary/read-processing387. The published supplement originals were attempted by the standard fetcher but remain behind PMC download gates in this pass. The mandatory scraper was run on the available article; displayed oligos were verified in its Methods.

This original implementation uses an expressed guide-barcode proxy (GBC), not direct capture of sgRNA. Antiparallel U6-sgRNA and EF1α-Puro/T2A/BFP/polyadenylated-GBC expression cassettes reside on the vector. A separately determined sgRNA↔GBC dictionary translates the RNA barcode to perturbation identity. The cell barcode (CBC), UMI, GBC and sequencer sample index have different roles.

## Source-positive local graph

🟢 The study loads sorted cells into the historical10xChromium system and follows its original manufacturer instructions. Main Methods describes Drop-seq testing as additional unpublished work, not the library chemistry of the selected10x study. The local graph begins with source whole-transcriptome amplified cDNA; it does not reconstruct that historical frontend from a current kit.

🟢 Approximately5ng of WTA material enters25cycle dial-out PCR with Q5mastermix. Forward primer columns concatenate as P5 +8base sequencing barcode +printed R1 handle +BFP-targeting arm `TAGCAAACTGGGGCACAAGC`. The selected printed index member is `TCGCCTTA`; the model imports canonical P5/Nextera motifs but not a family recipe. The source reverse primer is `CAAGCAGAAGACGGCATACGAGAT` (P7). Complete forward transcription is checked against the actual source text. Exact landing on the original historical cDNA remains unassembled.

🟢 The source program gives98°C10s followed by25cycles containing98°C2s,65°C5s,72°C10s and final72°C1min. Product is run on2%agarose and the appropriate band extracted; NextSeq sequencing follows. Source describes paired-end reads containing CBC/UMIs on one side and GBC on the other, but the full polymer/read-primer geometry and run-cycle allocation are not yet assembled. The graph therefore stops at gel-purified material and claims no acquired sequencing endpoint.

## Separate dictionary library

The paper also prepares a paired-end MiSeq library from intermediate plasmids to associate sgRNA and GBC before cassette insertion. Those primer/index tables are not used as if they were the post-WTA GBC capture library. Dictionary assignment and duplicate/ambient-barcode filtering are computational operations, not molecular product changes. GBC width and full vector sequence are unassigned here rather than invented from another Perturb-seq implementation.

## Deferrals and component ancestry

Final-VOR comparison, published chemistry/vector supplements, full GBC transcript architecture, historical10x guide and exact final read layout remain open. No missing-information claim or independent absence pass is asserted. Explicit source ancestry candidates: historical10x capture backend (Methods318–322), and the vector shared with Adamson et al.2016 (Methods256); component/version scope must be checked before importing either recipe. Direct-capture Perturb-seq and CROP-seq remain distinct protocols.

## Completion

Bounded original-source pass only. Catalogue remains notes. Private partial graphs retain explicit uncovered source inventory and unresolved terminal handoffs; passing structural/local checks does not complete source review. No two-reader missing-information pass has been performed.

## Bounded publication review, 30 October 2026

The associated accepted manuscript Methods318–405 was re-read, including the original primer table, the actual 50 µl Q5 reaction, gel recovery, NextSeq statement and paired-end processing. Mandatory whole-bundle oligo scraping was repeated; the source table was checked directly. This is one original, not a published-VOR comparison or an absence audit.

The current model remains a private post-WTA enrichment fragment. It has no acquired sequencing terminal and retains the required `unassembled_source` item. Its existing source-positive 5 ng WTA allocation, indexed dial-out PCR and 2% agarose recovery do not themselves account for that declared remaining input/library scope. The original historical capture/RT preparation, final occupied strands and primer/read geometry remain unassembled. This is a model/source-scope hold, not a claim that the paper omitted these details. The selected WTA-to-GBC boundary may be useful independently, but its exact input identity and declared source scope still need reconciliation; dictionary preparation or a modern Chromium recipe must not fill the gap by analogy.

The actual article delegates the cell-to-WTA preparation to historical manufacturer instructions. The relevant document/revision is not identified in this pass; no unavailable disclosure document is invented for an unassigned component. The manifest retains real gated supplementary originals. Existing source requests for FigureS1/FigureS5/TableS1 and the final-VOR comparison remain open. No fresh retrieval or source expansion occurred.
