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
