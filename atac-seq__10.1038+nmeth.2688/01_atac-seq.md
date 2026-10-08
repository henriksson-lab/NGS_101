# ATAC-seq — research notes

## 1. What it is

ATAC-seq profiles accessible chromatin by allowing adapter-loaded Tn5 to cleave and tag
native chromatin in one reaction.

## 2. Sources read

- 🟢 Buenrostro et al. 2013, DOI [10.1038/nmeth.2688](https://doi.org/10.1038/nmeth.2688),
  full author manuscript and Online Methods (PMC3959825).
- 🟢 The paper's Supplementary Table 1 for Ad1_noMX, Ad2.1_TA, Nextera Ad_noMX and indexed
  PCR primer design.

## 3. Construct-changing steps

🟢 Mild lysis releases nuclei. Nextera Tn5 is incubated with native chromatin, directly
fragmenting accessible DNA and transferring its mosaic-end adapters. 🟢 A 72 °C extension
precedes five PCR cycles; a qPCR side reaction determines the additional cycle count.

## 4. Final library and reads

🟡 The productive two-ended molecule is the standard heterologous Nextera product and is
therefore built from the shared canonical components. Paired reads enter the genomic insert
from opposite transposition sites; i5 and i7 are PCR-added sample indexes.
