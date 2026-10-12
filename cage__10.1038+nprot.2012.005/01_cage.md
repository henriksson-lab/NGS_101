# CAGE

## 1. What it is

CAGE enriches capped RNA 5′ ends by cap trapping and converts them into short sequencing tags for promoter and transcription-start-site mapping.

## 2. Authoritative source

🟢 Takahashi et al., *5′ end-centered expression profiling using cap-analysis gene expression (CAGE) and next-generation sequencing*, Nature Protocols (2012), DOI `10.1038/nprot.2012.005`.

## 3. Construct-changing operations

🟢 First-strand cDNA remains paired to RNA while the cap is oxidized and biotinylated. RNase digestion and streptavidin capture select cap-complete hybrids before linker installation and tag construction.

## 4. Exact oligos and read geometry

🟢 Table 1 gives RT-N15-EcoP as `AAGGTCTATCAGCAGNNNNNNNNNNNNNNN`, PCR forward as `AATGATACGGCGACCACCGACAGGTTCAGAGTTC`, PCR reverse as `CAAGCAGAAGACGGCATACGA`, and the custom sequencing primer as `CGGCGACCACCGACAGGTTCAGAGTTCTACAG`.

🟢 The 5′ linker contains a three-base sample barcode, `CAGCAG`, and a terminal N6/GN5 ligating tract. The 3′ linker begins with a two-base protrusion. EcoP15I releases a 27-nt cap-derived tag.

🟡 The final duplex is assembled from the overlapping PCR-primer and linker sequences. The custom run primer is located against that same construct rather than assumed to be standard TruSeq.

## Ligation geometry

🟢 For the AGA representative, the N6 upper strand is `CCACCGACAGGTTCAGAGTTCTACAGAGACAGCAGNNNNNN` with a 3′ phosphate; the GN5 alternative replaces the first random position with G. The matching lower is `/5Phos/CTGCTGTCTCTGTAGAACTCTGAACCTGTCGGTGG/3AmMO/`. The fixed 35-base portions form an exact antiparallel duplex. The upper strand's terminal N6/GN5 tail anneals transiently at a captured cDNA end, aligning the lower strand's 5′ phosphate beside the cDNA 3′-OH. This is a source-defined splinted ligation, not an unknown linker end.

🟢 The 3′ upper linker is `NNTCGTATGCCGTCTTCTGCTTG` and the lower is `CAAGCAGAAGACGGCATACGA`. The initial `NN` is therefore a **5′ two-base protrusion** on the upper strand. The earlier schematic placed it at the 3′ end; the oligo table fixes both its position and polarity.

🟢 The defining runs were single-read (36 cycles on GAIIx or 50 cycles on HiSeq 2000). The custom sequencing primer first reads the inline sample barcode and EcoP15I site, then the cap-derived transcript tag.
