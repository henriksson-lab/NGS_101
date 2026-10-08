# NEBNext Single Cell / Low Input RNA

## 1. What it is

A commercial full-length cDNA workflow for individual cells or very low RNA inputs, followed by Ultra II FS Illumina library conversion.

## 2. Sources and evidence

- 🟢 NEB E6420 manual and [product page](https://www.neb.com/products/e6420-nebnext-single-cell-low-input-rna-library-prep-kit-for-illumina).
- 🟢 The product uses an RT primer, template-switching oligo and cDNA PCR primer to generate amplified full-length cDNA.
- 🟢 Amplified cDNA is processed by Ultra II FS fragmentation/end preparation, adaptor ligation and indexed PCR.
- 🔴 NEB names but does not publish the RT primer, TSO or cDNA PCR-primer sequences; the schematic preserves those boundaries as unknown.

## 3. Molecular path

Template switching supplies the second amplification handle without ligating an adaptor to the RNA. The resulting full-length cDNA is amplified first; only then is it fragmented for an ordinary paired-end sequencing library.
