# ARTIC / PrimalSeq tiled-amplicon sequencing

## 1. Protocol boundary

This is the two-pool viral tiling protocol defined by Quick et al. in 2017, before later
SARS-CoV-2-specific ARTIC schemes and kit revisions. A primer scheme is a replaceable input; the
protocol's transferable molecular idea is assigning neighboring overlapping amplicons to alternate
multiplex reactions and then attaching platform-specific libraries.

## 2. Source read

🟢 Quick et al., *Multiplex PCR method for MinION and Illumina sequencing of Zika and other virus
genomes directly from clinical samples*, *Nature Protocols* (2017),
DOI `10.1038/nprot.2017.066`.

The complete protocol PDF listed in `ref/MANIFEST.md` was read, including Figure 3 and steps 6–12.
The published Supplementary Information (`nprot.2017.066_41596_2017_BFnprot2017066_MOESM1_ESM.pdf`)
was also checked: Tables 1 and 2 print the full ZikaAsian and ChikAsiaECSA
400-nt primer schemes, respectively.

## 3. Tiled multiplex PCR

🟢 Primal Scheme chooses overlapping amplicons across a viral reference. Alternate regions go to
pool 1 and pool 2, so neighboring overlapping primers never share a reaction; otherwise a short
overlap product can outcompete the intended products.

🟢 Each primer is used at 0.015 µM. PCR uses 98 °C denaturation and a five-minute 65 °C combined
annealing/extension step. Each pool is bead-cleaned separately. For the paper's approximately
400-base design, paired 250-base MiSeq reads can cover each product; MinION accepts much longer
amplicons.

🟢 The first ZikaAsian pool-1 pair is `ZIKA_400_1_LEFT`
`GACAGTTCGAGTTTGAAGCGAAAG` and `ZIKA_400_1_RIGHT`
`AGTATGCACTCCCACGTCTAGT` (Supplementary Table 1). These are target-specific
PCR primers, not universal library adapters. The supplement supplies the remaining
ZikaAsian and ChikAsiaECSA pairs with their pool assignments; the page shows one
representative pair and keeps the scheme replaceable.

## 4. MinION branch

🟢 Products are normalized by molecular amount, end-repaired and dA-tailed. Native barcodes are
ligated, barcoded samples are pooled, and the sequencing adapter is ligated to the pool. The paper
used the 2D EXP-NBD002/SQK-LSK208 workflow and also states compatibility with the then-current 1D
EXP-NBD103/SQK-LSK108 workflow. Because those kits changed frequently, the authors deliberately
did not specify their internal component names or oligo sequences.

## 5. MiSeq branch

🟢 Each pool contributes 50 ng. KAPA Hyper performs end repair and dA-tailing, then SureSelectXT2
indexing adapters replace the KAPA adapters during ligation. The protocol recommends paired
250-base sequencing for 400-base amplicons.

🟡 The paper does not print the SureSelectXT2 oligo bases. The schematic marks a canonical
Illumina-compatible final adapter geometry as inferred while keeping the reported paired-end and
sample-index read layout explicit.
