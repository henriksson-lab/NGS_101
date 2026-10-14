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
4. 🟢 Prepare both pulldown and input libraries with KAPA HyperPrep reactions: blunt
   end, 3′ dA tail, ligate adaptors carrying six-base indexes, and amplify 8–16 cycles
   with KAPA HiFi polymerase.

## Library and run

🟢 The paper names KAPA HyperPrep and six-base indexed adaptors but does not identify
the adaptor supplier or print their oligo sequences. The contemporary
[KAPA HyperPrep technical sheet](https://www.encodeproject.org/documents/1bca5eae-b9ce-4df0-b401-612924dc53a4/@@download/attachment/03_KAPA_Hyper_Prep_Kit_KR0961-TDSv2.15_final_24%20rxn_Reddy_GGR.pdf)
explicitly states that adapters are not included in the kit. Therefore the paper's
KAPA kit name cannot establish the indexed adaptor identity. The compatible
Illumina shell on the schematic remains inferred.

🟢 Libraries were pooled and run for 101 cycles on HiSeq 2000, 2500 or 4000 instruments,
yielding 100-nt single-end reads. The reported TSA-seq score is the normalized pulldown
read density relative to normalized input read density across genomic windows.
