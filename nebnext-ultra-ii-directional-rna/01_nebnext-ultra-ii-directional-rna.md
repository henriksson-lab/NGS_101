# NEBNext Ultra II Directional RNA

## 1. What it is

A commercial strand-specific bulk RNA-seq workflow based on dUTP marking and selective removal of the second cDNA strand.

## 2. Sources and evidence

- 🟢 [NEB E7760/E7765 manual v7.0 (April 2025)](https://www.neb.com/en-us/-/media/nebus/files/manuals/manuale7760_e7765.pdf), read from `/husky/henriksson/NGS101_papers/sources/nebnext-ultra-ii-directional-rna/manuale7760_e7765.pdf.txt`, and the [official RNA workflow](https://www.neb.com/en-us/products/next-generation-sequencing-library-preparation/library-preparation-for-illumina/rna-library-prep-for-illumina).
- 🟢 NEB states that dUTP is incorporated during second-strand synthesis and that USER excises this labelled strand.
- 🟢 The workflow additionally supplies Strand Specificity Reagent.
- 🟢 The E7760/E7765 kit does **not** include the adaptor or amplification primers: its manual requires a separate NEBNext Oligo Kit or customer-supplied adaptor and primers (manual page 2). A specific oligo kit must therefore be named before exact indexed primers can be assigned to this workflow.
- 🟢 The manual prints Read 1 `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` and Read 2 `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` as **adapter trimming sequences** (page 3). These are sequences encountered in reads, not an oligo ordering sheet; they do not specify the ligation-ready adapter's complete strands or modifications.
- 🟡 Exact random-primer identities and kit-internal sequence boundaries are not disclosed; the page models their roles without assigning bases.

## 3. Molecular path

Fragmented RNA is random-primed, copied to first-strand cDNA, converted to a dUTP-bearing duplex, end-prepared and adapter-ligated. USER destroys the marked second strand before indexed PCR, so reads retain transcript orientation.
