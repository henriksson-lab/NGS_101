# NEBNext Ultra II Directional RNA

## 1. What it is

A commercial strand-specific bulk RNA-seq workflow based on dUTP marking and selective removal of the second cDNA strand.

## 2. Sources and evidence

- 🟢 NEB product E7760/E7765 manual and [official RNA workflow](https://www.neb.com/en-us/products/next-generation-sequencing-library-preparation/library-preparation-for-illumina/rna-library-prep-for-illumina).
- 🟢 NEB states that dUTP is incorporated during second-strand synthesis and that USER excises this labelled strand.
- 🟢 The workflow additionally supplies Strand Specificity Reagent.
- 🟡 Exact random-primer identities and kit-internal sequence boundaries are not disclosed; the page models their roles without assigning bases.

## 3. Molecular path

Fragmented RNA is random-primed, copied to first-strand cDNA, converted to a dUTP-bearing duplex, end-prepared and adapter-ligated. USER destroys the marked second strand before indexed PCR, so reads retain transcript orientation.
