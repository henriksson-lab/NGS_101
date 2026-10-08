# NEBNext Enzymatic Methyl-seq v2

## 1. What it is

An enzymatic alternative to bisulfite sequencing that reports 5mC and 5hmC together while avoiding bisulfite-mediated DNA damage.

## 2. Sources and evidence

- 🟢 NEB E8015 manual v2.0, June 2025, and [product page](https://www.neb.com/en-us/products/e8015-nebnext-enzymatic-methyl-seq-v2-kit).
- 🟢 TET2 oxidizes 5mC through 5caC; T4-BGT converts 5hmC to glucosylated 5hmC.
- 🟢 APOBEC deaminates unmodified C to U but does not convert the protected products.
- 🟢 The manual prints both TruSeq-compatible adaptor trimming sequences and specifies Q5U amplification after conversion.

## 3. Molecular path

Fragmented DNA is end-repaired, dA-tailed and adapter-ligated before protection and deamination. After Q5U amplification, genomic unmodified C is read as T, whereas protected 5mC and 5hmC are read as C; the assay does not distinguish those two modified states.
