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

## 3. Molecular order

🟢 Cells undergo one round of replication in BrdU. Each chromatid therefore contains one
unsubstituted parental strand and one BrdU-containing nascent strand. A single daughter
cell is isolated. Hoechst sensitization and UV introduce nicks selectively into BrdU DNA;
the damaged nascent strand is excluded during library construction. Sequencing reads from
the retained templates map in chromosome-scale orientation patterns.

## 4. Historical library architecture

🟢 The Online Methods specify Illumina PE adapters, PE 1.0, a custom P7/index/PE2 PCR
primer, paired 76-nt reads, and the custom index-read primer sequence. The barcode is
called a fault-tolerant hexamer even though the displayed degenerate run has seven Ns;
the seven-cycle index read therefore contains the six barcode bases plus the adjacent
fixed P7′ base. 🟡 The standard PE adapter and Read 1/Read 2 primer bases are matched to
Illumina's authoritative obsolete paired-end oligo table, rather than to modern TruSeq.
