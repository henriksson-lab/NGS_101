# NEBNext Single Cell / Low Input RNA

## 1. What it is

A commercial full-length cDNA workflow for individual cells or very low RNA inputs, followed by Ultra II FS Illumina library conversion.

## 2. Sources and evidence

- 🟢 NEB E6420 manual and [product page](https://www.neb.com/products/e6420-nebnext-single-cell-low-input-rna-library-prep-kit-for-illumina).
- 🟢 The E6420 manual's "Oligo Sequences" table (p. 23, Rev 6.0) prints all three cDNA-stage reagents, including the three terminal riboguanosines on the TSO:
  - RT primer mix: `5′-AAGCAGTGGTATCAACGCAGAGTAC(T)30V-3′` (`V` is A/C/G). The table writes out all 30 T bases.
  - Template switching oligo: `5′-GCTAATCATTGCAAGCAGTGGTATCAACGCAGAGTACATrGrGrG-3′`.
  - cDNA PCR primer: `5′-AAGCAGTGGTATCAACGCAGAGT-3′`.
- 🟢 Amplified cDNA is processed by Ultra II FS fragmentation/end preparation, adaptor ligation and indexed PCR.
- 🟢 The same manual prints the Read 1 and Read 2 adaptor-trimming sequences, `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` and `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`.
- 🔴 The exact kit adaptor body, sample-index sequences and corresponding library PCR primers are not in that table. The trimming sequences are partial library motifs, not the complete adapter oligos.

## 3. Molecular path

Template switching supplies the second amplification handle without ligating an adaptor to the RNA. The resulting full-length cDNA is amplified first; only then is it fragmented for an ordinary paired-end sequencing library.
