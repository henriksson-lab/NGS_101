# SureSelect XT HS2 DNA target enrichment

## 1. What it is

Agilent SureSelect XT HS2 converts fragmented genomic DNA into dual-indexed Illumina libraries, optionally places a five-base molecular barcode at each insert end, pools indexed libraries before capture, hybridizes them to biotinylated RNA baits, retrieves bait–library duplexes on streptavidin beads and performs post-capture PCR.

## 2. Authoritative source

🟢 Agilent, *SureSelect XT HS2 DNA Kits: Library Preparation (+/−MBC), Overnight-Hyb Target Enrichment, Pre-capture Pooling Workflow*, G9957-90000 version B0, May 2026. The manual explicitly specifies mechanical or enzymatic fragmentation, combined end repair/dA-tailing, adaptor ligation, pre-capture indexing PCR, blocker-plus-probe hybridization, streptavidin capture, post-capture PCR and standard Illumina paired-end sequencing.

🟢 The manual states that MBC-tagged libraries carry inline molecular barcodes at both insert ends and instructs analysis software to trim five bases from the beginning of each read.

## 3. Construct logic

1. Fragment genomic DNA mechanically or enzymatically.
2. Repair ends and add 3′ dA overhangs.
3. Ligate the HS2 adaptor mix; the MBC version places a five-base tag next to each insert end.
4. PCR adds the paired eight-base P5/P7 indexes and completes the Illumina library.
5. Pool indexed libraries, block their adaptor arms and hybridize target intervals to complementary biotinylated RNA probes overnight.
6. Streptavidin beads retain probe-bound molecules; stringent washes remove non-target molecules.
7. Post-capture PCR replenishes the selected library without changing its read structure.

## 4. Uncertainty and scope

🟡 The current manual publishes index sequences but not the complete adaptor-oligo sequences. The schematic therefore uses canonical Illumina sequencing arms and marks the kit-specific adaptor duplex as a role rather than inventing hidden bases. Capture baits are design-specific and are shown as a complementary target interval.

## 5. Read layout

Read 1 and Read 2 each encounter a five-base MBC before target DNA. Index 1 and Index 2 read the paired eight-base sample indexes. Primer placement is computed on the final construct.
