# NEBNext Enzymatic Methyl-seq v2

## 1. What it is

An enzymatic alternative to bisulfite sequencing that reports 5mC and 5hmC together while avoiding bisulfite-mediated DNA damage.

## 2. Sources and evidence

- 🟢 [NEB E8015 manual v2.0, July 2026](https://www.neb.com/en/-/media/nebus/files/manuals/manuale8015.pdf), and [product page](https://www.neb.com/en-us/products/e8015-nebnext-enzymatic-methyl-seq-v2-kit). The local copy is `manuale8015.pdf.txt`, read in the 2026-10-17 bounded pass. Sections read: kit contents, Required Materials (p. 2), Overview and trimming sequences (pp. 2–3), protocol §§1.1–1.11 (pp. 4–12), Kit Components and Revision History (pp. 13–14). The revision history lists 1.2 (9/25), 1.3 (4/26) and 2.0 (7/26). The "June 2025" revision named in the initial note does not match any of these, and earlier revisions have not been compared.
- 🟢 TET2 oxidizes 5mC through 5caC; T4-BGT converts 5hmC to glucosylated 5hmC.
- 🟢 APOBEC deaminates unmodified C to U but does not convert the protected products.
- 🟢 The manual says the libraries "resemble TruSeq libraries" and prints two trimming sequences (p. 3). Adaptor Read 1 is `AGATCGGAAGAGCACACGTCTGAACTCCAGTCA` and Adaptor Read 2 is `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT`. The self-test checks that both read-through sequences occur in the drawn library.
- 🟢 Workflow (§§1.1–1.10):
  - Combine 48 µl sample DNA (0.1–200 ng) with 1 µl each of diluted unmethylated lambda and CpG-methylated pUC19 controls.
  - Fragment to ~350 bp with NEBNext UltraShear or a Covaris, with no clean-up.
  - Run Ultra II End Prep at 20 °C for 15 min, then 65 °C for 15 min.
  - Ligate the NEBNext EM-seq Adaptor with Ultra II Ligation Master Mix and Ligation Enhancer, 20 °C for 15 min. Clean up with 1X beads; at ≤10 ng input, add Carrier DNA.
  - Protect with TET2, T4-BGT (diluted 1:10 at ≤10 ng), UDP-glucose and DTT. Diluted Fe(II) starts the 1 h reaction at 37 °C. Add Stop Reagent and incubate 30 min at 37 °C. Clean up with 1X beads.
  - Denature with formamide (recommended) or 0.05 N NaOH, 85 °C for 10 min, then a chilled block.
  - Deaminate with APOBEC: 37/42/50/75 °C steps at ≤100 ng, or 37 °C for 3 h for any input. Go straight to PCR without clean-up.
  - Amplify with Q5U Master Mix and an LV UDI primer pair for 4–14 cycles (98/62/65 °C).
  - Clean up with 0.8X beads.
  - Sequence at a user-chosen length (typically 2×76, 2×100 or 2×150).
- 🟡 Model corrections in this pass:
  - Added the denaturation step, which the model had omitted.
  - Added the control spike-ins, fragmentation, Stop Reagent and the three bead clean-ups.
  - The adaptor-ligated intermediate no longer carries P5/P7 and index windows. Those are added by the UDI PCR.
  - The adaptor end is now declared unavailable. The manual does not print a 3′-dT end or call the adaptor "conversion-resistant"; earlier text asserted both.
- 🟢 The E8015 manual specifies separately supplied NEBNext LV Unique Dual Index Primer Sets (E3390/E3392/E3400/E3402/E3404/E3406/E3408). [Their own manual](https://www.neb.com/en/-/media/nebus/files/manuals/manuale3390_e3392_e3400_e3402_e3404_e3406_e3408) publishes the expected i7 and i5 read sequences by well. Thus the sample-index bases can be assigned for a chosen plate/well; the E8015 assay itself has no single fixed index pair.
- 🟢 **Published discrepancy:** E8015 v2.0_7/26 p. 2, “Required Materials Not Included,” names UDI kits **E3390/E3392**. The Fig. 1 caption on the same page prints **E3900/E3902**. The dedicated LV UDI manual title and index tables also name **E3390/E3392**; the model follows those kit identifiers. NEB has not confirmed a correction to the Fig. 1 caption, and the June 2025 revision has not been compared.
- 🔴 The E8015 and linked LV UDI manuals inspected so far print index-read values by well. They do not disclose the complete EM-seq Adaptor strands, the adaptor's ligating end or cytosine chemistry, or complete PCR-primer synthesis oligos. The trimming sequences contain C and are read as C after conversion, which implies the adaptor cytosines are not deaminated. The manual does not say how (🟡). This missing-information review remains open pending an independent check of the newly recovered UDI manual and the earlier E8015 revision.

## 3. Molecular path

Fragmented DNA is end-repaired, dA-tailed and adapter-ligated before protection and deamination. After Q5U amplification, genomic unmodified C is read as T, whereas protected 5mC and 5hmC are read as C; the assay does not distinguish those two modified states.
