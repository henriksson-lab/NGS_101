# TSA-seq — research notes

Chen et al., *Journal of Cell Biology* (2018),
[doi:10.1083/jcb.201807108](https://doi.org/10.1083/jcb.201807108).

## Source read

- 🟢 Defining article and Materials and methods (PMCID: PMC6219710; JATS text twin
  Results 86–102, Materials and methods 226–262: TSA cell labeling, DNA isolation,
  fragmentation, dot blot, affinity purification, library preparation and sequencing).
- 🟢 KAPA Hyper Prep Kit technical data sheet KR0961 v2.15 (ENCODE-hosted copy), product
  description lines 136–165 and Table 2 (ER & AT reaction). The paper does not cite a
  sheet revision, so this copy describes the kit, not the exact lot used.

## What the assay records

🟢 Antibody-localized horseradish peroxidase creates short-lived tyramide-biotin radicals.
Their steady-state concentration decays away from the labelled nuclear compartment, and
DNA reacts directly with those radicals. Consequently the probability that a genomic
fragment carries biotin is a population-level function of cytological distance. There is
no proximity-ligation junction, spatial barcode or per-molecule distance value.

## Construct-changing steps

1. 🟢 Formaldehyde-fix (10 ml 8% PFA added to 40 ml culture, 20 min RT, glycine quench) and permeabilize
   (0.5% Triton X-100) cells, quench endogenous peroxidases with 1.5% H2O2 for 1 h, block,
   stain the chosen compartment marker with primary antibody (20–24 h, 4 °C) and an
   HRP-conjugated secondary (20–24 h, 4 °C), then perform tyramide-biotin signal
   amplification for 10 min (8 mM tyramide-biotin stock diluted 1:5,000 with 0.003% H2O2;
   conditions 2 and 3 add 50% sucrose, and condition 3 also DTT). A no-primary control was
   processed in parallel; it is not drawn as a separate branch. The defining demonstrations targeted SON-marked speckles and
   nuclear lamins.
2. 🟢 Lyse with SDS/proteinase K at 55 °C, reverse crosslinks in 0.2 M NaCl at 65 °C for
   20–24 h, phenol–chloroform extract, RNase A treat and ethanol precipitate DNA. Sonicate to 200–800 bp, with most fragments
   around 500 bp, and reserve 10% as input.
3. 🟢 Bind biotinylated DNA to Dynabeads M-270 streptavidin, wash stringently (W&B, TSE
   500 and TE washes; a later protocol used 55 °C W&B/Triton washes), digest the beads
   overnight with proteinase K in the presence of free biotin, phenol–chloroform extract
   and precipitate with glycogen carrier. The paper estimates only ~0.002–0.15% of
   ~500-bp fragments carry biotin, so the unbound fraction is mostly unlabelled DNA.
4. 🟢 Prepare both pulldown and input libraries with the Kapa Hyper Prep kit: blunt
   end and 3′ dA tail, ligate indexed adaptors with 6-nt barcodes, and amplify 8–16 cycles
   with KAPA HiFi polymerase. 🟡 The kit sheet performs end repair and A-tailing in one
   ER & AT reaction, so the model draws one combined step. Post-ligation cleanups are
   not described by the paper and are not drawn.

## Library and run

🟢 The paper names KAPA HyperPrep and six-base indexed adaptors but does not identify
the adaptor supplier or print their oligo sequences. The contemporary
[KAPA HyperPrep technical sheet](https://www.encodeproject.org/documents/1bca5eae-b9ce-4df0-b401-612924dc53a4/@@download/attachment/03_KAPA_Hyper_Prep_Kit_KR0961-TDSv2.15_final_24%20rxn_Reddy_GGR.pdf)
explicitly states that adapters are not included in the kit. Therefore the paper's
KAPA kit name cannot establish the indexed adaptor identity. The compatible
Illumina shell on the schematic remains inferred.

🟢 Libraries were pooled and run for 101 cycles on HiSeq 2000, 2500 or 4000 instruments,
yielding 27–81 million 100-nt reads per sample. 🟡 Single-end acquisition and a 6-cycle
index read are model inferences from the 101-cycle count and the 6-nt barcodes; the
paper does not state the read mode or index-read cycles. The reported TSA-seq score is the normalized pulldown
read density relative to normalized input read density across genomic windows.

## Open questions

- Adaptor supplier and sequences, amplification primers and any post-ligation cleanup
  used by the UIUC sequencing core (not identified in the inspected sources).
- Actual read mode and index-read cycles.
- No-primary control branch and staining conditions 1–3 are not drawn separately.
