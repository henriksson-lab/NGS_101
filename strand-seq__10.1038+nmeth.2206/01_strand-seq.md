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

## 4. Evidence boundary

🟢 Strand selection and its inheritance signal are defining-paper facts. 🟡 Historical
library-kit bases are unavailable in the article, so the page marks canonical paired-end
arms as inferred while keeping the actual selected molecule explicit.
