# Ribo-ITP

Ozadam et al.2023, DOI https://doi.org/10.1038/s41586-023-06228-9 . Private capped source-qualified review, no finished publication status or accepted absence claim.

🟢 Published PMC10307641 main overview70–90, Methods188–248, processing256; original oligo workbook MOESM9 read. The reporting-summary PDF MOESM2 selected identity/data-collection section was also read; it is not a detailed chemistry supplement. Those are three originals; other named supplements and the authors' online protocol are deferred at the cap. Mandatory scraper ran over the fetched bundle, skipping its large data workbook; actual adopted bases/modifications were verified separately in the oligo original.

🟢 Microfluidic isotachophoresis extracts and size-selects ribosome-protected RNA from nuclease-digested lysate; it does not create RNA sequence by itself. Fluorescent marker DNA has actual terminal ddC nucleotide plus a3′ extension block, not a removable text annotation. Source D-Plex small-RNA library preparation is modified with PNK and diluted TSO; half the cDNA is amplified, leaving an independent remainder. RNA-seq comparator and conventional monosome-isolation preparations remain separate.

## Physical preparation and library

🟢 The published main supplies the actual preparation despite deferred detailed online protocol/buffer supplements. Individual K562 cells are sorted using singlet/DAPI-negative gates into cold cycloheximide lysis buffer and flash frozen. The single-cell MNase branch digests37°C30min and quenches with EGTA; the RNaseI branch digests37°C15min and quenches with SDS. The source100-cell and mouse-oocyte/embryo preparations have their own printed collection/dilution details, retained as modes rather than extra methods. Conventional bulk sucrose-gradient monosome isolation and comparator Smart-seq3/SMARTer recipes are not inherited.

🟢 Fluorescent markers are added to digested lysate before electrophoretic extraction. Nucleic acids travel through5% and10% polyacrylamide stages; the elution reservoir is washed and exchanged to dephosphorylation buffer. The between-marker10µl fraction is collected. Other channel/reservoir/wash material remains a physical companion; no perfect marker removal, purity or size cutoff is predicted. The separate demonstration that collected earlier/later fractions is not conflated with the default retained fraction.

🟢 D-Plex small-RNA kitC05030001 is modified with0.5µl T4PNK and25min dephosphorylation. Fig2c explicitly describes poly(A) tailing and template-switching RT; Methods dilute TSO1:2. Kit-internal RT/primer/polymer geometry stays unreviewed. Half of cDNA is amplified17cycles, leaving a separate mixture remainder with unresolved later fate; no storage instruction is inferred. A1:4 library dilution is analytically assessed with BioAnalyzer; assayed material does not feed the production pool. Source approximate equimolar pooling explicitly includes other molecular library inputs (membership/sequences unassigned),1.8xAMPure and source size-selection branches follow before NovaSeq6000 sequencing. Run cycles, complete end/primer geometry and exact index members remain unassigned.

🟢 All100-cell human samples and three MII oocytes use the single-index module; other mouse samples use unique dual indexes. Source BluePippin windows differ accordingly. The human single-cell index-module assignment is not filled in by analogy. RNaseI single-cell libraries use PAGE purification of200bp products. Source UMI extraction reads12nt at the5′ read end and discards four following RT-derived nucleotides computationally. This is not a physical deletion from the library; the UMI region has structured identity metadata, with full source oligo placement unassembled.

## Modified oligo source

Published oligo workbook MOESM9 is named TableS7 by the main article; its own header says Supplementary Table5. Source identity is anchored by the exact manifest URL and oligo rows, not a guessed table number. Original XLSX shared-string values were inspected independently of extraction. Active DNA marker strings as written:

```
/5ATTO488N/CTGCTGGAGTTCGTGACC/3ddC/
/5ATTO488N/CCTCCTTGTATAAATCCTGGTTGCTGTCTCTTTAT/3ddC/
```

The terminal `/3ddC/` contributes a cytidine nucleotide with dideoxy sugar; it is retained in the actual marker sequence and as an extension block. The source's functional statement is that markers do not amplify in downstream library preparation. The local constructor checks only the narrower selected reaction requiring3′OH; it does not claim every conceivable polymerase or reaction cannot ever process the molecule. Dephosphorylating sampleRNA does not silently turn a DNA ddC into a normalOH terminus.

Other workbook RNA controls (17nt ddC-ended,26nt and34nt3′phosphorylated) were inventoried as original-written modified RNA; they are not active DNA marker substitutions or universal RPF sequences. The source guard matches the represented active DNA marker polymers and whole modified transcription is tested locally. Synthetic `SelectedRPF` expectations establish only sequence retention/dephosphorylation/tailing relations, not a full kit recipe.

## Limits

Read originals: published main, selected reporting-summary identity/data-collection section and oligo workbook. The reporting summary is not a chemistry supplement. Detailed buffer Table6/MOESM8, other figure/analysis supplements, the authors' online Ribo-ITP protocol and historical D-Plex module guides remain deferred at the three-original cap. Large MOESM7 data workbook was not scanned by the mandatory scraper; this is recorded, not treated as absent sequence evidence. Full template/end/index/read geometry, actual carryover and exact downstream fates remain unresolved. No private schematic, public promotion or accepted absence claim.

Independent selected-source implementation review corrected two provenance assumptions: other DNA libraries are explicit molecular inputs to pooling, and the unused cDNA half has UNRESOLVED disposition rather than unsupported storage. Focused checks cover both.

## Bounded unlinked-page adjudication (2026-10-30, batch24)

EGTA and SDS are now physical molecular inputs with separate quench conditions. The source RT supplies yield functional catalyst, nucleotide, primer/TSO and buffer inventories; PCR supplies similarly yield actual functional components. Packaging, number of source stocks, doses and whole historical kit oligos remain unassembled. These are provenance decompositions, not a premix factory or modern D-Plex recipe. The existing half-cDNA remainder and analytical Bioanalyzer aliquot stay separate from production. Actual retained library/pool/cleanup products now retain library typing; this does not prove whole sequence assembly or perfect purity. The source 10 µl selected eluate is collected from dephosphorylation buffer, without an invented extra 3 µl addition.

The original `pending` item remains required and uncovered: online protocol, buffer Table6, relevant SI, historical D-Plex inventory/assembly and read geometry are still unreviewed. Published main and original MOESM9 scope stay separate from historical Reporting Summary identity-only scope. No public promotion or missing-information claim.
