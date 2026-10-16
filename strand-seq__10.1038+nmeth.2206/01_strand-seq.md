# Strand-seq — research notes

## 1. What it is

Strand-seq selectively sequences the parental template DNA strands inherited by one
daughter cell, retaining chromosome-wide Watson/Crick orientation information.

## 2. Sources read

- 🟢 Falconer et al. 2012, DOI
  [10.1038/nmeth.2206](https://doi.org/10.1038/nmeth.2206), complete author manuscript
  (PMC3580294).
- 🟢 Sanders et al. 2017 / later protocol descriptions for the explicit Hoechst/UV
  BrdU-strand destruction mechanism; used as clarification rather than as defining source.
- 🟢 The defining paper's archived Supplementary Figures and Supplementary Data files were checked for an oligo table; they show library and chromosome results but do not supply the individual fault-tolerant hexamer barcode set.

## 3. Molecular order

🟢 Cells undergo one round of replication in BrdU. Each chromatid therefore contains one
unsubstituted parental strand and one BrdU-containing nascent strand. A single daughter
cell is isolated. Hoechst sensitization and UV introduce nicks selectively into BrdU DNA;
the damaged nascent strand is excluded during library construction. Sequencing reads from
the retained templates map in chromosome-scale orientation patterns.

## 4. Historical library architecture

🟢 The Online Methods specify fragment end repair, A tailing, Illumina PE adapter
ligation, PE 1.0, a custom P7/index/PE2 PCR primer, paired 76-nt reads, and the custom
index-read primer sequence. In “DNA template strand library construction” the barcode
is called a fault-tolerant **hexamer** although the printed PCR primer has **seven Ns**;
“Illumina sequencing” specifies a **seven-base** index read. 🟡 The model uses six variable
index bases plus one adjacent fixed base, but the paper does not explicitly resolve the
seven-N printing. The standard PE adapter and Read 1/Read 2 primer bases are matched to
[Illumina's obsolete paired-end oligo table](https://support-docs.illumina.com/SHARE/AdapterSequences/Content/OligoSeqsPairedEndDNA.htm), rather than to modern TruSeq.
The PE adapter's 3′ dT end pairs with the source-described A-tailed fragment. 🔴 The
particular hexamer assigned to each of the 66 libraries has not been verified from the
linked supplementary material.
