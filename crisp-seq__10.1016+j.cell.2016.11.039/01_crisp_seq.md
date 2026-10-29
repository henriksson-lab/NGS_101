# Dissecting Immune Circuits by Linking CRISPR-Pooled Screens with Single-Cell RNA-Seq

Defining DOI: https://doi.org/10.1016/j.cell.2016.11.039

## Defining original and distinct preparation boundary

🟢 Jaitin etal., Cell167:1883–1896.e1–e7(15December2016), DOI10.1016/j.cell.2016.11.039. Institutional CELL2016.pdf printed title/DOI/journal pagination verified. Read STARMethods vector constructionpp e2–e3 (text2510–2540), CRISP library preparationp e4 (2603–2640) and UGI processing2657–2663. Mandatory scraper rerun after download; complete library oligos were checked in the original text. No co-publication or application-only argument is used for admission.

The molecular distinction is source-specific: a polyadenylated fluorophore transcript contains an8nt Unique Guide Identifier(UGI) afterWPRE and372bp upstream of the SV40polyA signal. The UGI↔gRNA association is established by Sanger sequencing. The guideRNA itself is not captured as a polyadenylated transcript. More importantly, the UGI readout is produced from an UNfragmented aliquot of MARS-seq amplifiedRNA with a specific adapter/secondRT/extra pool-barcodePCR. This differs from original Perturb-seq's post-WTA DNA dial-out PCR, CROP-seq's LTR architecture and direct guideRNA capture.

## Upstream identity and physical branch split

🟢 Sorted-cell mRNA is simultaneously barcoded/converted to cDNA, pooled and T7-amplified under the cited historical MARS-seq procedure. DNaseI treatment and1.2×SPRI recovery precede the local graph input. That source frontend is not reconstructed from an unspecified modern kit.

🟢 Half of the aRNA enters a fragmented MARS-seq transcriptome-library branch. Separately10–12% enters UGI-seq WITHOUT fragmentation. Remaining material has unresolved fate, not automatic disposal. The local graph preserves these distinct physical branches and stops the gene branch after fragmentation because its ligation/RT/PCR recipe remains unassembled.

## UGI-specific oligos and prefix

🟢 UGI ligation adaptor: `ATGATCAAGCGACCACCGAG`,5′phosphate and3′C3 spacer. The C3 blocker is non-nucleotide; no extra cytidine is appended. Original RNA-end processing and ligase formulation are inherited by citation only, not read/validated here, so complete ligation geometry remains unavailable.

🟢 After source cleanup, secondRT uses `CTCGGTGGTCGCTTGATCAT`, the computed full reverse complement of the adapter body. The consumed adapter reagent is one strand carrying the source5′phosphate and3′C3 blocker in its actual model state. A separate binding view adds the later RT primer; that primer is not silently present in the ligation reagent. Computed pairing is not full RNA ligation validation. Local unresolved-geometry drawings say unreviewed/unassembled, not unreported by the source. No undocumented RNA hydrolysis is inserted.

🟢 Intermediate10cyclePCR adds a pool/plate barcode using:

`CTACACGACGCTCTTCCGATCTNNNNNXXXXTCCCCGCGTCGACGGATC`

where sourceN5 means random bases with unassigned function and sourceXXXX is a4base platebarcode. N5 is NOT called a new UMI. Reverse P7_Rd2 is `CAAGCAGAAGACGGCATACGAGATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`. Final MARS-likePCR remains unassembled; source provides NextSeq500 readout but no model-acquired run geometry.

## Identifier roles, deferrals and ancestry

The source processing distinguishes cellbarcode7, original RMT/UMI8 and UGI8, with4base platebarcode introduced later. These source roles do not fix unassembled full-strand positions or transform the intermediateN5 window into a UMI. UGI decoding, quality thresholds and molecule collapsing are computational; no rejected read is modeled as physically removed DNA.

Exact historical MARS-seq end processing/ligation, original RT/cell/UMI template, intermediate workup, finalPCR/index/read layout and linked detailed originals remain deferred. The paper explicitly cites Jaitin2014 MARS-seq for the frontend/library operations, with a local unfragmented UGI branch and source-specific oligos. This is component-level evidence for ancestry; no source-selected shared constructor is imported and no absence claim is made.
