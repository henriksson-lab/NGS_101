# CHANGE-seq — research notes

Lazzarotto et al., *Nature Biotechnology* (2020),
[doi:10.1038/s41587-020-0555-7](https://doi.org/10.1038/s41587-020-0555-7).

## Assay target

The sequenced material is circularized human genomic DNA linearized in vitro at Cas9 cleavage sites. This protocol does not target, amplify or reconstruct a pathogen genome: every oligo it uses is a host-genome or library oligo (uracil-bearing custom Tn5 ends and library primers). The note is here so that the protocol's CRISPR, off-target and double-strand-break vocabulary is not read as pathogen work.

## Sources read

- 🟢 The defining article (PMC7652380 JATS, CHANGE-seq Methods paragraph, line 172) and its
  detailed Supplementary Protocol (MOESM1 text: reagent setup lines 100–130, procedure
  steps 1–30 lines 155–500, sequencing lines 575–589).
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
2. 🟢 Tn5 fragments purified genomic DNA to approximately 400 bp (55 °C, 7 min; step 2),
   proteinase K stops the reaction, and gap repair with 2X Kapa HiFi HotStart Uracil+ Ready
   Mix plus Taq DNA ligase at 72 °C for 30 min completes the staggered transposition
   products (steps 2–8, lines 170–211).
3. 🟢 One combined USER/T4 PNK reaction (37 °C, 1 h; step 9) exposes and phosphorylates the
   engineered ends, after which T4 DNA ligase at 16 °C for 16 h circularizes each fragment
   (step 11, lines 212–250); the article gives the circularization concentration as 5 ng/µl.
4. 🟢 Exonuclease I, lambda exonuclease and Plasmid-Safe ATP-dependent DNase (37 °C 1 h,
   70 °C 30 min; step 13) destroy residual molecules with free ends. Closed circles survive
   (lines 262–290).
5. 🟢 Re-folded sgRNA and SpCas9 (article: 90 nM Cas9, 270 nM sgRNA on 125 ng circularized
   DNA, 37 °C 1 h; steps 15–16) cut only circles containing an on- or off-target site.
   Proteinase K then removes the RNP (step 17), leaving new linear DNA ends (lines 291–350).
6. 🟢 The selected products are dA-tailed with the Kapa A-tailing enzyme (30 °C, 30 min;
   step 19) and ligated to the NEBNext Adapter for Illumina by Kapa DNA ligase (20 °C, 1 h;
   step 21) — the hairpin adapter is NEB's, the tailing and ligation enzymes come from the
   Kapa HTP Library Preparation Kit PCR-free. USER (supplied with the NEBNext Dual Index
   Primers Set 1 kit) opens the hairpin (step 23), and NEBNext i5/i7 primers add dual indexes
   during 20 Kapa HiFi cycles (steps 25–27, lines 351–490).
7. 🟢 The Supplementary Protocol pools the libraries and sequences them on a MiSeq
   600-cycle v3 kit with 150-bp paired-end and 8-bp dual-index reads (step 40, line 589);
   the article's Methods paragraph reports MiSeq or NextSeq 550, 150-bp paired-end. Both
   genomic reads start at a Cas9-generated end and point into the selected circle-derived
   fragment (defining article Methods; Supplementary Protocol reagent list and workflow).

## What is and is not an identifier

The four-base custom-transposon palindrome is a cohesive circularization end, not a UMI
or sample barcode. CHANGE-seq uses the NEBNext i5 and i7 sample indexes. The defining
method does not assign a molecule-counting UMI.

## Open questions

- 🔴 The commercial NEBNext Set 1 index pair actually used for each library is not printed in
  the article or its supplements; the i501/i701 members shown on the page are vendor examples.
- 🟡 The exact end chemistry left by USER on the custom transposon (a four-base 3′ ACGT
  overhang derived from oCRL225/oCRL226) is read out of the two printed oligos rather than
  drawn in the source.

## Interpretation boundary

CHANGE-seq measures biochemical cleavage of purified genomic DNA. Circularization and
exonuclease treatment enrich newly Cas9-linearized molecules extremely efficiently, but
the assay does not establish that a candidate site is accessible or edited inside a cell.

## Source-owned validation migration, 2026-10-18

Local graph now owns the E7600S Set1 endpoint, eight-cycle dual indexes, Kapa PCR inputs and exact binding geometry on the USER-opened adapter template. The template is the preceding open_hairpin product, not the final indexed library. Canonical vendor arms/primer architecture remain kit-qualified, with index members unassigned. Custom dU transposon consumes a paired phosphorylated oligo scene; dA-tail state uses explicit opposed strands/overhangs instead of two dA labels on one strand. Distinct symbolic intermediates preserve opened ends/circle/cut provenance without claiming exact sequence derivation. Tagmentation off-product recovery remains unresolved instead of an unsupported physical discard.

Read source scope: Published Supplementary Information 41587_2020_555_MOESM1_ESM.pdf.txt: reagent/oligo section lines57–105, library steps19–27 lines351–480 and sequencing steps39–40 lines585–590; published article PMC7652380 Methods, extracted line172.

Remaining work: Actual study-specific Set1 member assignment and physical kit reagent verification beyond the named product; Exact input-output relations for Tn5 gap repair, USER opening, circle closure and nuclease linearization remain unvalidated; graph cleanups remain incomplete. No accepted missing-information claim is added.
