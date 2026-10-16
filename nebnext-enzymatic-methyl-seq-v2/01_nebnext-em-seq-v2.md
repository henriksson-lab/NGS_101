# NEBNext Enzymatic Methyl-seq v2

## 1. What it is

An enzymatic alternative to bisulfite sequencing that reports 5mC and 5hmC together while avoiding bisulfite-mediated DNA damage.

## 2. Sources and evidence

- 🟢 [NEB E8015 manual v2.0, July 2026](https://www.neb.com/en/-/media/nebus/files/manuals/manuale8015.pdf), and [product page](https://www.neb.com/en-us/products/e8015-nebnext-enzymatic-methyl-seq-v2-kit). The earlier June 2025 revision cited in the initial note remains unavailable for direct comparison.
- 🟢 TET2 oxidizes 5mC through 5caC; T4-BGT converts 5hmC to glucosylated 5hmC.
- 🟢 APOBEC deaminates unmodified C to U but does not convert the protected products.
- 🟢 The manual prints both TruSeq-compatible adaptor trimming sequences and specifies Q5U amplification after conversion.
- 🟢 The E8015 manual specifies separately supplied NEBNext LV Unique Dual Index Primer Sets (E3390/E3392/E3400/E3402/E3404/E3406/E3408). [Their own manual](https://www.neb.com/en/-/media/nebus/files/manuals/manuale3390_e3392_e3400_e3402_e3404_e3406_e3408) publishes the expected i7 and i5 read sequences by well. Thus the sample-index bases can be assigned for a chosen plate/well; the E8015 assay itself has no single fixed index pair.
- 🟢 **Published discrepancy:** E8015 v2.0_7/26 p. 2, “Required Materials Not Included,” names UDI kits **E3390/E3392**. The Fig. 1 caption on the same page prints **E3900/E3902**. The dedicated LV UDI manual title and index tables also name **E3390/E3392**; the model follows those kit identifiers. NEB has not confirmed a correction to the Fig. 1 caption, and the June 2025 revision has not been compared.
- 🔴 The E8015 and linked LV UDI manuals inspected so far print index-read values by well but do not disclose complete conversion-resistant adaptor strands or complete PCR-primer synthesis oligos. This missing-information review remains open pending an independent check of the newly recovered UDI manual and the earlier E8015 revision.

## 3. Molecular path

Fragmented DNA is end-repaired, dA-tailed and adapter-ligated before protection and deamination. After Q5U amplification, genomic unmodified C is read as T, whereas protected 5mC and 5hmC are read as C; the assay does not distinguish those two modified states.
