# LamPORE

## Defining source

🟢 James *et al.*, “LamPORE: rapid, accurate and highly scalable molecular screening for SARS-CoV-2 infection, based on nanopore sequencing”, medRxiv (2020), [doi:10.1101/2020.08.07.20161737](https://doi.org/10.1101/2020.08.07.20161737).

The later clinical implementation paper, [doi:10.1128/JCM.03271-20](https://doi.org/10.1128/JCM.03271-20), confirms the 96-well combinatorial layout.

## Molecular workflow

🟢 Four multiplexed RT-LAMP primer sets target SARS-CoV-2 ORF1a, E and N plus human ACTB. One of eight sample barcodes is embedded in each forward inner primer (FIP), so it is copied repeatedly into the branched LAMP products.

🟢 Eight differently FIP-barcoded reactions are pooled by plate column. A rapid transposase library step adds one of twelve second barcodes, yielding 8 × 12 combinatorial sample identities. The motor-adapted fragments are read directly by Nanopore sequencing.

🔴 The commercial FIP-barcode, rapid-barcode, attachment-tag and motor-adapter sequences are not disclosed by the defining sources. The page therefore represents their molecular roles and does not invent bases or one linear consensus for heterogeneous LAMP concatemers.

## Read interpretation

🟢 Classification uses the rapid barcode, the repeatedly copied FIP barcode, and target-derived sequence assigning a read to ORF1a, E, N, ACTB or unclassified. Nanopore acquisition is primer-free and variable-length.
