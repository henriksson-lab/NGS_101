# NicE-seq — research notes

Ponnaluri et al., *Genome Biology* (2017),
[doi:10.1186/s13059-017-1247-6](https://doi.org/10.1186/s13059-017-1247-6).

## Sources read

- 🟢 The defining open-access article and methods (PMCID: PMC5488340).
- 🟢 The authors' GEO record GSE97889, particularly GSM2579632/SRX2741932, for
  deposited-library metadata and paired-end run identity.
- 🟢 New England Biolabs, *NEBNext Ultra DNA Library Prep Kit for Illumina*, E7370
  manual, for the discontinued commercial end-prep and ligation workflow named by the
  paper. The paper does not identify which separately sold NEBNext oligo set was used.

## What carries the accessibility signal

1. 🟢 Cells are either used natively or crosslinked in 1% formaldehyde for 10 min and
   quenched with 125 mM glycine. NP-40 cytosolic buffer isolates nuclei.
2. 🟢 Nt.CviPII recognizes `CCD` (`D` = A, G or T) and nicks one strand where chromatin
   permits access. The new 3′-OH primes *E. coli* DNA polymerase I. Polymerase extension
   plus 5′→3′ exonuclease displacement translates the nick through a tract of DNA.
3. 🟢 The 200-µl labeling reaction contains 2.5 U Nt.CviPII, DNA polymerase I and 30 µM
   each dNTP, including 6 µM biotin-14-dATP and 6 µM biotin-16-dCTP. It runs for 2 h at
   37 °C with 800-r.p.m. mixing. Biotin is therefore internal to accessible DNA, not an
   adapter-end label.
4. 🟢 EDTA stops the reaction while RNase A removes RNA. Proteinase K and SDS reverse
   crosslinks at 65 °C, and phenol/chloroform extraction recovers covalently labelled
   genomic DNA.
5. 🟢 Covaris sonication makes approximately 150-bp fragments. The E7370 workflow
   end-repairs, dA-tails and ligates an Illumina-compatible adapter *before* enrichment.
6. 🟢 Adapter-ligated DNA is incubated with blocked streptavidin beads for 2 h in 2 M
   NaCl. Four B&W/Triton washes and a TE/Triton wash remove unlabelled molecules. PCR
   from the bead-bound fraction, usually ten cycles, yields the sequencing library.

This order matters: read endpoints are mechanical shear positions and do not identify
the original Nt.CviPII nick. Enrichment across the fragment reports proximity to an
internally biotinylated accessible-chromatin tract.

## Source discrepancy

🟢 The article Methods state 10 U DNA polymerase I in the standard labeling reaction,
whereas the authors' GSM2579632 extraction metadata state 50 U for that deposited sample.
The public schematic deliberately says “DNA polymerase I” without resolving this conflict
into a false universal amount. Both sources agree on 2.5 U Nt.CviPII, nucleotide and
biotin-nucleotide concentrations, temperature, time and mixing.

## Commercial-oligo boundary

🟢 E7370 supplies end-repair/dA/ligation/PCR reagents but requires a separately purchased
NEBNext oligo kit. The article says only “NEB Illumina adaptor (E7370S)” and does not name
the oligo catalogue number, sample index or whether a later dual-index option was used.

🟡 The schematic consequently shows a compatible historical six-base, single-i7
TruSeq-shaped NEBNext shell as **inferred**. It does not claim an exact index sequence and
does not invent an i5 read. The biochemical heart of NicE-seq—accessible-site nicking,
biotin nick translation, random shearing and internal-label capture—is unaffected by that
unresolved purchasing choice.

## Sequencing

🟢 The deposited NicE-seq data are paired-end NextSeq 500 libraries and were processed as
read pairs. The article does not state a NicE-seq-specific cycle count (its explicit 75-bp
statement is for the comparison ATAC-seq library), so no R1/R2 length is assigned on the
public page. Standard TruSeq-family primers bind a compatible E7370 library shell.
