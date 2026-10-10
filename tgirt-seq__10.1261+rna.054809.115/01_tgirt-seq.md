# TGIRT-seq — source notes

Defining application: Qin et al., *RNA* 2016, [doi:10.1261/rna.054809.115](https://doi.org/10.1261/rna.054809.115). Detailed current protocol and oligo table: Xu et al., *Bio-protocol* 2021, [doi:10.21769/BioProtoc.4239](https://doi.org/10.21769/BioProtoc.4239).

## Sources read

- 🟢 PMC8678547 HTML text, Figure 1 description, procedure, and complete oligonucleotide list.
- 🟢 PMC4691826 full text for the defining plasma-RNA implementation and its original molecular layout.

## Evidence trail

- 🟢 A 35-nt, 3′-blocked R2 RNA anneals to a 36-nt R2R DNA primer carrying a mixed one-base 3′ overhang.
- 🟢 TGIRT-III template switches through the single base pair between the overhang and target-RNA 3′ nucleotide, thereby joining the R2R sequence to the cDNA while reverse transcription proceeds.
- 🟢 Alkali degrades RNA; the completed cDNA is purified before a 5′-adenylated, 3′-blocked R1R DNA is attached by thermostable 5′ App DNA/RNA ligase.
- 🟢 An optional six-base random UMI precedes R1R. PCR primers add P5/P7 and optional i5 plus required i7 indices; the protocol recommends dual indices for NovaSeq.
- 🟢 The current R2 RNA/R2R pair contains the documented extra A/T near the cDNA end to reduce R1R–R2R dimers.
- 🟢 The final page deliberately selects the printed 6N-UMI R1R option and the protocol's six-base TSBC01-style barcode example on both PCR primers. Other barcode choices and omission of i5 or UMI are also permitted by the protocol but are not silently merged into this construct.

## Oligos as ordered

- R2 RNA: `5′-rArArG rArUrC rGrGrA rArGrA rGrCrA rCrArC rGrUrC rUrGrA rArCrU rCrCrA rGrUrC rArC/3SpC3/-3′`
- R2R DNA: `5′-GTG ACT GGA GTT CAG ACG TGT GCT CTT CCG ATC TTN-3′`
- 6N UMI R1R: `5′-/5Phos/NNN NNN GAT CGT CGG ACT GTA GAA CTC TGA ACG TGT AG/3SpC3/-3′`
- P5 PCR primer: `5′-AAT GAT ACG GCG ACC ACC GAG AT BARCODE C TAC ACG TTC AGA GTT CTA CAG TCC GAC GAT C-3′`
- P7 PCR primer: `5′-CAA GCA GAA GAC GGC ATA CGA GAT BARCODE GTG ACT GGA GTT CAG ACG TGT GCT CTT CCG ATC T-3′`
