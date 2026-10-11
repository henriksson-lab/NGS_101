# CHANGE-seq — research notes

Lazzarotto et al., *Nature Biotechnology* (2020),
[doi:10.1038/s41587-020-0555-7](https://doi.org/10.1038/s41587-020-0555-7).

## Sources read

- 🟢 The defining article and its detailed Supplementary Protocol, including the reagent
  setup, exact custom transposon oligos and steps 2–27.
- 🟢 The NEBNext Multiplex Oligos for Illumina (Dual Index Primers Set 1) manual already
  transcribed into `lib/illumina.py`; the defining protocol names kit E7600S but does not
  reprint its commercial oligos.
- 🟢 Supplementary Tables 2–8 were inspected to distinguish guide/validation primers from
  oligos that become part of the CHANGE-seq assay library.

## Construct logic

1. 🟢 The custom transposome is annealed from oCRL225,
   `5′-/5Phos/ACG/ideoxyU/AGATGTGTATAAGAGACAG-3′`, and oCRL226,
   `5′-/5Phos/CTGTCTCTTATACACATCTACGT-3′` (Supplementary Protocol lines 101–119).
   The paper describes a 19-bp Tn5 mosaic end plus a four-base palindromic end.
2. 🟢 Tn5 fragments purified genomic DNA to approximately 400 bp. Uracil-tolerant gap
   repair and Taq ligase complete the staggered transposition products (lines 156–211).
3. 🟢 USER and T4 PNK expose and phosphorylate the engineered ends, after which dilute
   overnight T4-ligase treatment circularizes each fragment (lines 212–250).
4. 🟢 Exonuclease I, lambda exonuclease and Plasmid-Safe ATP-dependent DNase destroy
   residual molecules with free ends. Closed circles survive (lines 262–290).
5. 🟢 Cas9–sgRNA cuts only circles containing an accessible on- or off-target site.
   Proteinase K removes the RNP, leaving new linear DNA ends (lines 291–350).
6. 🟢 The selected products are dA-tailed and ligated to the NEBNext Adapter for Illumina.
   USER opens the hairpin adapter, and NEBNext Set 1 i5/i7 primers add dual indexes during
   20 PCR cycles (lines 351–490).
7. 🟢 Libraries are sequenced with 150-bp paired-end reads on MiSeq or NextSeq. Both
   genomic reads start at a Cas9-generated end and point into the selected circle-derived
   fragment (defining article Methods; Supplementary Protocol reagent list and workflow).

## What is and is not an identifier

The four-base custom-transposon palindrome is a cohesive circularization end, not a UMI
or sample barcode. CHANGE-seq uses the NEBNext i5 and i7 sample indexes. The defining
method does not assign a molecule-counting UMI.

## Interpretation boundary

CHANGE-seq measures biochemical cleavage of purified genomic DNA. Circularization and
exonuclease treatment enrich newly Cas9-linearized molecules extremely efficiently, but
the assay does not establish that a candidate site is accessible or edited inside a cell.
