# ARTIC / PrimalSeq tiled-amplicon sequencing

## Assay target

Unlike the other break-mapping and screen-readout protocols documented here, this one does sequence a viral genome: tiled amplicons across a clinical Zika or chikungunya sample. It is a diagnostic and genomic-surveillance method, and its primer schemes are published in full by the authors and maintained openly by the ARTIC Network for that purpose. What this page documents is the library chemistry downstream of the two primer pools; the per-amplicon primer coordinates are a replaceable published scheme file, not anything derived here.

## 1. Protocol boundary

This is the two-pool viral tiling protocol defined by Quick et al. in 2017, before later
SARS-CoV-2-specific ARTIC schemes and kit revisions. A primer scheme is a replaceable input; the
protocol's transferable molecular idea is assigning neighboring overlapping amplicons to alternate
multiplex reactions and then attaching platform-specific libraries.

## 2. Source read

🟢 Quick et al., *Multiplex PCR method for MinION and Illumina sequencing of Zika and other virus
genomes directly from clinical samples*, *Nature Protocols* (2017),
DOI `10.1038/nprot.2017.066`.

The published typeset article (*Nat Protoc* 12:1261–1276, co-author mirror listed in
`ref/MANIFEST.md`; local `quick-multiplex-zika-sequencing.pdf.txt`) is the authority. The 2026-10-17
bounded pass read Experimental design (primer design, multiplex PCR, library preparation;
pp. 1263–1266), Figure 3 legend (p. 1266), Materials (p. 1267) and Procedure Steps 4–12 (pp. 1268–1272).
The PMC copy (PMC5902022) is the NIHMS author manuscript; it is superseded and not used for
active claims. The published Supplementary Information (`nprot.2017.066_41596_2017_BFnprot2017066_MOESM1_ESM.pdf`)
was also checked: Tables 1 and 2 print the full ZikaAsian and ChikAsiaECSA
400-nt primer schemes, respectively.

## 3. Tiled multiplex PCR

🟢 Primal Scheme chooses overlapping amplicons across a viral reference. Alternate regions go to
pool 1 and pool 2, so neighboring overlapping primers never share a reaction; otherwise a short
overlap product can outcompete the intended products.

🟢 Q5 Hot Start High-Fidelity polymerase is used with each primer at 0.015 µM (Step 6 table and
critical step; 36 primers in ZikaAsian pool 1). An earlier repository note reported a 0.02 µM
prose value. That value is in the author manuscript only, not in the published article, so it is a
manuscript-to-publication revision and not a discrepancy. Step 8 cycling is 98 °C for 30 s, then
cycles 2–40 of 98 °C for 15 s and 65 °C for 5 min. The published primer-design section describes the
two-step protocol with 95 °C denaturation, which differs from the cycling table. The page shows this
as a published discrepancy, and the model uses 98 °C. Each pool is cleaned separately with AMPure XP
(1.0×, 0.8× or 0.6× by amplicon length; Step 9) and kept separate. For the paper's approximately
400-base design, paired 250-base MiSeq reads can cover each product; MinION accepts much longer
amplicons.

🟢 The first ZikaAsian pool-1 pair is `ZIKA_400_1_LEFT`
`GACAGTTCGAGTTTGAAGCGAAAG` and `ZIKA_400_1_RIGHT`
`AGTATGCACTCCCACGTCTAGT` (Supplementary Table 1). These are target-specific
PCR primers, not universal library adapters. The supplement supplies the remaining
ZikaAsian and ChikAsiaECSA pairs with their pool assignments; the page shows one
representative pair and keeps the scheme replaceable.

## 4. MinION branch

🟢 Products are normalized by molecular amount (Step 12A(ii) table, ~0.3 pM total input per flow cell)
and end-repaired/dA-tailed with NEBNext Ultra II End Prep (20 °C for 5 min, then 65 °C for 5 min).
Native barcodes NB01–NB12 are ligated with Blunt/TA Ligase Master Mix (20 °C for 10 min, then 65 °C
for 10 min). One barcode per pool per sample is recommended, and one barcode per sample is allowed.
The barcoded reactions are pooled and SPRI-cleaned, and the barcoding adapter is ligated per the kit
protocol, with Blunt/TA ligase in place of Quick Ligation. The graph draws the recommended
per-pool barcoding. The paper
used the 2D EXP-NBD002/SQK-LSK208 workflow and also states compatibility with the then-current 1D
EXP-NBD103/SQK-LSK108 workflow. Because those kits changed frequently, the authors deliberately
did not specify their internal component names or oligo sequences.

## 5. MiSeq branch

🟢 Each pool is kept in its own tube (50 ng), and two barcodes are recommended per sample. KAPA Hyper
performs end repair and dA-tailing, and SureSelectXT2 indexing adapters (G9622A) replace the KAPA
adapters during ligation. A 0.8× SPRI cleanup replaces the 1× post-amplification cleanup "to remove
potential adaptor–dimers", and the graph shows these as discarded candidate off-products. The
libraries then go through TapeStation QC and KAPA qPCR quantification and are pooled equimolar.
Sequencing uses MiSeq Reagent Kit v2 (500 cycle), and the protocol recommends paired 250-base reads
for 400-base amplicons.

🟡 The paper does not print the SureSelectXT2 oligo bases. The schematic marks a canonical
Illumina-compatible final adapter geometry as inferred while keeping the reported paired-end and
sample-index read layout explicit.

## 6. Open questions (2026-10-17 bounded pass)

- The KAPA Hyper Prep (07962363001) and SureSelectXT2 index (G9622A) originals are unread. The
  adapter, ligation-end and library-PCR primer bases are therefore not modeled.
- The ONT EXP-NBD103/SQK-LSK108 and EXP-NBD002/SQK-LSK208 kit protocols are unread.
- The ARTIC protocol site versions are unread.
- The 95 °C versus 98 °C denaturation statements are both in the published article and are not
  reconciled.
