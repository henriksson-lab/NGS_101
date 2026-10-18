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
   amplification for 10 min (buffer B contains 8 mM tyramide-biotin stock diluted 1:5,000 with 0.003% H2O2, then is mixed with an equal volume of buffer A;
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
   with KAPA HiFi polymerase. The native graph distinguishes the source-named end-repair, dA-tailing, ligation and PCR transitions; it assigns no additional post-ligation cleanup. The historical kit-sheet reading is not a fresh review of the study lot.

## Library and run

🟢 The inspected Methods name KAPA HyperPrep and six-base indexed adaptors. Native adapter and amplification-primer identities remain unassigned in this bounded model. The contemporary
[KAPA HyperPrep technical sheet](https://www.encodeproject.org/documents/1bca5eae-b9ce-4df0-b401-612924dc53a4/@@download/attachment/03_KAPA_Hyper_Prep_Kit_KR0961-TDSv2.15_final_24%20rxn_Reddy_GGR.pdf)
explicitly states that adapters are not included in the kit. Therefore the paper's
KAPA kit name cannot establish the indexed adaptor identity. The previous compatible Illumina shell has been removed; a kit name does not establish that native sequence.

🟢 Libraries were pooled and run for 101 cycles on HiSeq 2000, 2500 or 4000 instruments,
yielding 27–81 million 100-nt reads per sample. Acquisition mode and index cycles remain unassigned: neither single-end acquisition nor a six-cycle index read follows from those values. The reported TSA-seq score is the normalized pulldown
read density relative to normalized input read density across genomic windows.

## Open questions

- Adaptor supplier and sequences, amplification primers and any post-ligation cleanup
  used by the UIUC sequencing core (not identified in the inspected sources).
- Actual read mode and index-read cycles.
- No-primary control branch and staining conditions 1–3 are not drawn separately.


## Bounded native migration (2026-10-18, batch07)

Inspected 1 originals: Published PMC6219710 Methods222–260: TSA labelling conditions, isolation, fragmentation, affinity-purification original/modified procedures, library preparation and sequencing. Mandatory bundle scraper reviewed; large source-data workbook skipped by scraper, not interpreted as absence.

Local graph replaces shared endprep/endpoint defaults. Source condition1 and later modified all-at-once capture are procedural references with actual sample pairing unassigned. DNA/protein/carrier populations are independent. Physical extraction/washes,10%input split, binding/release, ER/dA/ligation and actual adapted PCR inputs are explicit. DNA-biotin and original/carrier fates remain unresolved; unstained cells are not discarded. Six-base physical sample-index feature is separate from unassigned index acquisition;101cycles/100nt does not establish single-end mode. BufferB concentrations precede equal-volume A mixing.

Deferred: Native study-lot adapters, ligation geometry, amplification primers and exact acquisition mask; Sample-to-capture-variant mapping; no-primary and other staining conditions remain separate undeclared routes; Historical KAPA sheet was not freshly reread and does not establish the study adapter supplier. This capped pass does not establish complete review or source-wide absence.

## Independent implementation correction (2026-10-18)

Confirmed source biotin originals are not computational spatial barcodes, matched10percentinput and90percentpulldown remain distinct, actual native adapted substrates feedPCR, and physical6nt adapterbarcode is not asserted as6cycle acquisition. Source101cycles/100nt reads stay qualified without SE/PE overassignment. Corrected a definite carrier continuity gap: independent glycogen present before endprep now persists in the pulldown A-tailed Scene instead of disappearing and reappearing afterligation. Originalbiotin/carrier/offtarget fates stay unresolved.
