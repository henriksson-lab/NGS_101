# REAP-seq — bounded source notes

Peterson et al. (2017), *Multiplexed quantification of proteins and transcripts in single cells*, [published DOI](https://doi.org/10.1038/nbt.3973).

## Sources and scope

This pass uses the **published** supplementary discussion, figures and tables:
`nbt.3973_41587_2017_BFnbt3973_MOESM5_ESM.pdf` in the external source store.
Discussion lines 1–34; Figure 2 caption 194–203; Figure 3 caption 233–248;
Tables 6–7 at 3528–3678 were inspected. The rendered original Table 7 (PDF
page 35) was also read to resolve column assignment and read directions. Tables
and figures remain source originals, not inferred from scraper output.
The full-bundle oligo scraper was run and its relevant hits checked against these
locations. The main article institutional mirror returned HTTP 403 according to
the retrieval manifest; its complete Methods have **not** been read in this pass.
MOESM6 and the remainder of MOESM5 have not been comprehensively reviewed.
This is not an exhaustive missing-information audit or a complete protocol model.

## Published chemistry

1. 🟢 Cells are labeled with antibody-conjugated DNA. The antibody tag is a
   Nextera Read 1 sequence, an eight-base antibody identifier and a poly(dA)
   tail. Excess unbound tags are washed away before droplet processing.
2. 🟢 The source uses **10x Genomics scRNA-seq v1**, not an unspecified current
   Chromium kit. Its bead oligo contains a 14-base cell barcode, ten-base UMI
   and 30-base oligo(dT), with amplification/sequencing sequences.
3. 🟢 After lysis, poly(dT) primes both antibody DNA and cellular mRNA. The
   reverse transcriptase extends the DNA reporter through its DNA polymerase
   activity while generating mRNA-derived cDNA. The source describes a
   noncovalent antibody–dsDNA complex that persists until PCR.
4. 🟢 After breaking droplets and silane cleanup, 0.6× SPRI separates longer
   mRNA-derived cDNA (>~500 bp, bead fraction) from shorter antibody-derived
   dsDNA (~155 bp, supernatant). **The supernatant is retained** for the protein
   branch. It must not be globally discarded as ordinary scRNA-seq waste.
5. 🟢 Exonuclease I degrades free ssDNA in the protein fraction to reduce
   antibody/cell-barcode crosstalk. Protein dsDNA is amplified with P5/P7 and
   the source's sample-index primer. The RNA branch follows the cited 10x v1
   protocol; its complete chemistry is not reconstructed here.

## Oligos and manufacturing detail

🟢 Figure 2 prints the antibody tag's 5′ modification as `/5AmMC6/`.
Its Nextera sequence is `TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG`, imported from
`nextera.READ1_PRIMER` in the model. A concrete eight-base antibody barcode is
followed by 25 dA **except** that a barcode ending in A receives 24 dA to keep
consecutive A insertions below 26. The local `antibody_tag()` constructor enforces
that source-specific rule. The unknown barcode in the drawn example is not a
source-verified antibody identity.

🟢 Table 6 prints P5 as `AATGATACGGCGACCACCGAGATCTACAC` and P7 as
`ATCTCGTATGCCGTCTTCTGCTTG`. These are imported from `illumina.P5` and
`illumina.P7_RC`, respectively. The indexed P5 oligo is P5 + eight-base sample
index + `TCGTCGGCAGCGTCAGATGT` (20 bases of Nextera Read 1). The model uses the
published S517 value `GCGTAAGA` as a concrete index example, not a universal index.

The Table 6 P7 value is the reverse complement of the canonical forward P7 PCR
primer. Its intended strand use and amplification stages need reconciliation
against the main Methods. This pass preserves the source-written value and leaves
PCR geometry unresolved; it makes **no claim of a publication error**.

## Sequencing: antibody identity is acquired in i7

🟢 Table 7 defines the protein run:

| Run read | Cycles | Content / primer |
|---|---|---|
| Read 1 | 26 | Custom `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC`: cell barcode 14 bases, then 12 bases of P7; adapter sequence trimmed |
| i7 | 8 | Custom Nextera Read 1 primer reads the antibody barcode |
| i5 | 8 | Standard primer reads the sample index |
| Read 2 | 10 | Standard Read 2 reads the UMI |

🟢 A TruSeq Read 1 primer is also supplied in the custom Read 1 and i7 mixes to
sequence PhiX for diversity. It is not a second primer binding site on the protein
reporter and must not be forced to land on that construct.

🟡 The displayed 5′→3′ strand is reconstructed from Table 7's segment order and
read arrows: P5 → sample index → Nextera Read 1 → antibody barcode → poly(dA) →
UMI → reverse-complement TruSeq Read 2 → cell barcode → reverse-complement P7.
Every segment defaults to inferred in this partial model. The three represented
sequencing primers land computationally; the source-declared i5 read remains
unplaced because this pass has not resolved the instrument-specific standard
primer/orientation. This does not establish that an i5 sequence is undisclosed.

## Representation and deferred obligations

The partial graph starts with the protein-containing SPRI supernatant, not a cell.
It models Exonuclease I treatment and a PCR operation with explicit unavailable
geometry. An uncovered source item deliberately prevents publication-complete
validation: cell labeling/capture, dual DNA/RNA synthesis, fractionation and the
full RNA library branch remain obligations. The local preview draws only the
inferred protein strand and read layout; registration remains **notes**.

Deferred: main published article Methods; remaining supplement sections including
Table 1 antibody dictionary and bulk-assay Table 8; published MOESM6 identity and
relevance; cited Zheng et al. 10x v1 original and its specific kit documents;
source-written antibody conjugation and release chemistry; PCR stage/polarity
reconciliation; full RNA branch; platform-specific i5 orientation; independent
source reading. No sequence or reaction is declared confirmed absent.

## Lineage

🟢 Supplementary Discussion and Figure 3 explicitly identify 10x Genomics v1 as
the inherited RNA-processing and bead-primer platform and cite Zheng et al.,
*Massively parallel digital transcriptional profiling of single cells*,
10.1038/ncomms14049. This is component-scoped inheritance of the v1 capture/RNA
branch, not evidence for a current kit revision or CITE-seq ancestry. Similarity
to CITE-seq is a modality/architecture comparison only.

## 2026-10-30 bounded publication preflight

The current model remains a post-SPRI protein-branch reference. The normal publication gate still rejects its unaccounted `cell_to_libraries` source item: upstream antibody/cell labelling, droplet synthesis, the physical size split and full RNA branch are not represented. Eleven existing partial-reference checks passing do not satisfy that required coverage.

This bounded pass reinspected published MOESM5 Discussion/Tables 6–7 (original Table 7 PDF p35) and identified MOESM6 as the 23 August 2017 reporting summary, not a bench protocol. The main institutional PDF remains HTTP403 at `https://rnajc.ucsf.edu/sites/rnajc.ucsf.edu/files/nbt.3973_0.pdf`; no retrieval retry or modern 10x recipe substitution was attempted. Main Methods, applicable predecessor/vendor originals and exact P7/priming orientation remain deferred. The source item stays required, the existing model is unchanged and unpublished, and no information is claimed absent from publication.
