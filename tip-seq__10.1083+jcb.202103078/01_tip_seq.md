# TIP-seq

**Defining source:** doi:10.1083/jcb.202103078.

## Sources read

- 🟢 PMC full text and Methods (PMCID `PMC8600797`).
- 🟢 Supplementary Table S1, including the bulk and 384-member sciTIP ME-T7
  families, second-strand primers, sequencing primers and indexed PCR primers.

## Construct path

🟢 `ME-scT7-r5001` anneals to 5′-phosphorylated ME-rev only across the terminal
19-nt mosaic end. Its T7 promoter, r5 barcode and S5 handle remain a single-stranded
transposome overhang. pA–Tn5 transfers this end at antibody-bound chromatin and Taq
fills the staggered gap.

🟢 T7 RNA polymerase transcribes outward into the target-proximal genomic DNA.
Random-hexamer RT creates first-strand cDNA; the published `sss_sci-nXTv2` primer
then binds immediately downstream of the T7 start and creates the second strand.

🟢 A separate MeB-only Tn5 fragments the cDNA and supplies its opposite adapter.
Representative S502 and N701 primers then add the second-well i5/i7 indexes, P5 and
P7. The final molecule retains the first-well r5 barcode between the second-strand
site and S5.

🟡 The public page draws the representative r5001/S502/N701 combination. The source
table contains 384 r5 transposons and the full indexed-PCR families; the representative
does not imply that those are the only usable indexes.

🟢 Published Table S2, available as the extracted `TableS2.docx` supplement,
specifies the sciTIP run as Read 1 51 cycles, index 1 8 cycles, index 2 29
cycles and Read 2 51 cycles. It says the 29-cycle index 2 read carries r5,
the intervening spacer and i5, and gives bcl2fastq header handling for these
cycles. This confirms the model's read layout from the source, beyond the
oligo table alone.
