# TSA-seq — research notes

Chen et al., *Journal of Cell Biology* (2018),
[doi:10.1083/jcb.201807108](https://doi.org/10.1083/jcb.201807108).

## Source read

- 🟢 Defining article and Materials and methods (PMCID: PMC6219710).

## What the assay records

🟢 Antibody-localized horseradish peroxidase creates short-lived tyramide-biotin radicals.
Their steady-state concentration decays away from the labelled nuclear compartment, and
DNA reacts directly with those radicals. Consequently the probability that a genomic
fragment carries biotin is a population-level function of cytological distance. There is
no proximity-ligation junction, spatial barcode or per-molecule distance value.

## Construct-changing steps

1. 🟢 Formaldehyde-fix and permeabilize cells, stain the chosen compartment marker with
   primary antibody and an HRP-linked detection system, then perform tyramide-biotin
   signal amplification. The defining demonstrations targeted SON-marked speckles and
   nuclear lamins.
2. 🟢 Reverse crosslinks and isolate DNA. Sonicate to 200–800 bp, with most fragments
   around 500 bp, and reserve 10% as input.
3. 🟢 Bind biotinylated DNA to Dynabeads M-270 streptavidin, wash stringently, digest the
   beads with proteinase K in the presence of free biotin and recover DNA.
4. 🟢 Prepare both pulldown and input libraries with Kapa Hyper Prep: blunt end, 3′ dA
   tail, ligate adapters carrying six-base indexes, and amplify 8–16 cycles with the kit’s
   HiFi polymerase.

## Library and run

🟢 The paper names the commercial library kit and six-base indexed adapters but does not
publish their oligo sequences. The public schematic therefore marks its compatible
Illumina shell as inferred rather than presenting vendor-neutral adapter bases as a
transcription.

🟢 Libraries were pooled and run for 101 cycles on HiSeq 2000, 2500 or 4000 instruments,
yielding 100-nt single-end reads. The reported TSA-seq score is the normalized pulldown
read density relative to normalized input read density across genomic windows.
