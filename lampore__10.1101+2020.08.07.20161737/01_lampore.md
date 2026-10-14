# LamPORE

## Defining source

🟢 James *et al.*, “LamPORE: rapid, accurate and highly scalable molecular screening for SARS-CoV-2 infection, based on nanopore sequencing”, medRxiv (2020), [doi:10.1101/2020.08.07.20161737](https://doi.org/10.1101/2020.08.07.20161737).

The later clinical implementation paper, [doi:10.1128/JCM.03271-20](https://doi.org/10.1128/JCM.03271-20), confirms the 96-well combinatorial layout.

## Molecular workflow

🟢 Four multiplexed RT-LAMP primer sets target SARS-CoV-2 ORF1a, E and N plus human ACTB. One of eight sample barcodes is embedded in each forward inner primer (FIP), so it is copied repeatedly into the branched LAMP products.

🟢 Eight differently FIP-barcoded reactions are pooled by plate column. A rapid transposase library step adds one of twelve second barcodes, yielding 8 × 12 combinatorial sample identities. The motor-adapted fragments are read directly by Nanopore sequencing.

🔴 The commercial FIP-barcode, rapid-barcode, attachment-tag and motor-adapter sequences are not disclosed by the defining sources. The page therefore represents their molecular roles and does not invent bases or one linear consensus for heterogeneous LAMP concatemers.

🟢 The later clinical LamPORE paper (Peto *et al.*, PMC8315953, Methods) says its ORF1a, E and N target sets use previously published primers and cites Zhang *et al.*, [doi:10.2144/btn-2020-0078](https://doi.org/10.2144/btn-2020-0078). Table 1 of that paper's [bioRxiv version](https://doi.org/10.1101/2020.06.03.132894) prints six unbarcoded primers each for As1e (ORF1a), E1, N2 and ACTB. As examples, its F3 primers are As1e `CGGTGGACAAATTGTCAC`, E1 `TGAGTACGAACTTATGTACTCAT`, N2 `ACCAGGAACTAATCAGACAAG`, and ACTB `AGTACCCCATCGAGCACG`.

🟡 The cited set gives defensible target-binding **precursor** oligos, but the LamPORE product adds a sample barcode into each FIP and the clinical paper does not identify the complete resulting FIP molecules. The cited unbarcoded sequences should not be copied into a LamPORE ordering list as finished FIPs.

## Read interpretation

🟢 Classification uses the rapid barcode, the repeatedly copied FIP barcode, and target-derived sequence assigning a read to ORF1a, E, N, ACTB or unclassified. Nanopore acquisition is primer-free and variable-length.
