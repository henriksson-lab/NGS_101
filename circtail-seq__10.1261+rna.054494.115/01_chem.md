# circTAIL-seq — capped source-qualified notes

Gazestani et al. (2016), circTAIL-seq, a targeted method for deep analysis of RNA3′tails, reveals transcript-specific differences by multiple metrics. DOI10.1261/rna.054494.115.

🟢 Published main Results60–78/Methods150–160 and the published supplementary optimization Method/TableS1 were read. The original TableS1 on supplement p2 was rendered and checked for bold sample tags and degenerate-window widths. Mandatory bundle primer scraper completed; two originals inspected, no predecessor original followed.

This adapts circularRT-PCR to direct NGS amplicon preparation, without the cloning step of traditional individual-tail analysis. It circularizes RNA, not cDNA as in SEnd-seq, and uses gene-specific terminal-junction primers rather than genome-wide terminal capture. No chemistry is inherited merely from the similar circular topology.

🟢 One independent DNase-treated5µg totalRNA sample is circularized in200µl with T4RNA ligase/RNase inhibitor/25µMATP overnight at room temperature →phenol/chloroform pH5.2/overnight precipitation/wash/20µl →1µg recoveredRNA for20µl Episcript RT with gene-specific primers →separate diagnostic optimization and production allocations. Recovery yield is unquantified: the model does not infer that1µg equals4µl of the recovered sample or declare unused RNA chemically removed.

🟢 Published oligos, 5′→3′:

- CO1 cDNA primer: `GTAATGAGTACGTTGTAAAACTG`.
- CO3p cDNA primer: `AAACCCTCCTCCTCTGGCGGT`.
- Reverse Illumina primer: the source P5/read1 arm followed by `NNNNN` and `ACTTAAGCACACAAGACATAG` for CO1, or `GCCTTTTTATTCAATTTTAAACAATCC` for CO3p. A visual space inside the printed CO3p sequence is formatting, not a nucleotide.
- Forward primer: source P7, a six-base sample barcode, source read2 arm, diversity offset, and gene-specific3′end sequence. CO1 A/B barcodes are `CGATGT`/`TGACCA`, CO3p A/B `ACAGTG`/`CAGATC`. A uses `NNNN`, B `NNNNNN`. Target arms are `TTTTATATTTTTTTTTGACAAG` for CO1 and `AGGGTTTAGTTGGAATGAA` for CO3p.

Code imports canonical P5/P7/read arms, checks the whole source-written indexed oligos and attaches structured sample-index meaning to the6nt barcode. N4/N5/N6 are **diversity offsets**, as the original table states, not molecule UMIs or additional sample identifiers. A/B identify the source primer versions; the model does not infer which independently recovered biological specimen received each oligo.

🟢 Production uses KAPA2G Robust with weighted dNTP mixture; CO1 has33cycles and CO3p32cycles at95/62/72°C15s. The source says0.25µl cDNA per stated reaction but the scaling relation between optimization and100µl production is not resolved. The supplementary Method explicitly divides100µl into20µl aliquots and recombines after PCR; these are five parallel physical reactions in the graph, not five consecutive amplification rounds. The detailed diagnostic optimization uses independent aliquots and remains outside this selected production branch.

🟢 Products are selected by6%nativePAGE over the expected zero-tail to150nt-tail window, eluted using the cited Riley2014 method with400µl/two columns, ethanol precipitated and resuspended15µl. The source's external Riley component remains opaque. Separate Bioanalyzer and KapaQC fractions precede an approximately equal-quantity six-library pool. An independent PhiX library is then added for the underclustered MiSeqV2 paired-end150 run. The paper describes12barcoded samples/sixRNAs overall, while TableS1 supplies the two selected targets; other run-member identities remain unassigned.

🟡 Full nativeRNA sequences, intermolecular ligation alternatives, circularization/primer junction geometry, diagnostic allocation and run-member selection remain unassembled. Physical RNA/PCR remnants remain until an actual separation. Downstream trimming of diversity offsets, read merging and tail filtering never removes physical library molecules. Private incomplete graph only; no public promotion or accepted absence claim.

Deferred source scope:

- Riley2014 original recovery method, directly cited by main Methods152; modified400ulelution/twocolumns locally established, full component not inherited.
- Source-cited conventional circularRT-PCR predecessor originals, including Slomovic/Schuster; no original read within cap.
- Exact diagnostic optimization/production sample allocation, remaining run members and full transcript/junction geometries remain open questions.
