# TIP-seq

**Defining source:** doi:10.1083/jcb.202103078.

## Sources read

- 🟢 PMC full text and Methods (PMCID `PMC8600797`; JATS `.xml.txt` lines 165, 169,
  177 bulk TIP-seq, 187 sciTIP-seq, 195 demultiplexing).
- 🟢 Supplementary Table S2 (run setup and bcl2fastq handling).
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

## Published Index 2 discrepancy

🟢 Table S1 lists `nXTv2_i5_index_read` as
`AATGATACGGCGACCACCGAGATCTACAC` (the P5 sequence). Primed there, a forward-strand
index 2 read would give i5 at cycles 1–8 and r5 at cycles 24–29. 🟢 The Methods
sciTIP-seq note says instead that the index 2 read has the 6-bp r5 at positions 1–6
and the 8-bp i5 at positions 21–29, with a 15-bp spacer between them, on a NovaSeq
6000 SP v1.5 flow cell.

🟡 The model follows the Methods description of the acquired run. It uses a
reverse-complement-workflow Index 2 primer (the Nextera Index 2 sequence, reverse
complement of S5 + ME), which places r5 at cycles 1–6 and i5 at cycles 22–29. The protocol self-test
derives these spans from the declared library and run lengths.
The printed 21–29 span covers nine positions for an 8-bp index. The role of the Table S1
P5-sequence primer is still unexplained. The page shows this discrepancy.

## Reaction details from the bulk Methods (sciTIP uses the same steps)

🟢 The Methods give these steps in order:

1. Taq 5X Master Mix gap fill at 72 °C for 3 min, with the SPRI beads kept in the reaction.
2. HiScribe T7 IVT at 37 °C for 16–19 h, then 2.0× SPRI.
3. Random-hexamer priming and SMART MMLV RT.
4. RNase H at 37 °C for 20 min to degrade the RNA in RNA–cDNA hybrids.
5. Second-strand priming at 65 °C, then Taq at 72 °C for 8 min.
6. MeB-only Tn5 at 55 °C for 6 min, stopped with 4 M GuHCl.
7. NEBNext High-Fidelity indexed PCR, with cycles set by qPCR (typically 7–12 for sciTIP).
8. A final 0.85× SPRI left-side size selection (>200 bp).

The graph names these enzymes. It does not draw the SPRI cleanups, GuHCl stripping or
size selection.

## Open questions

- Bulk TIP-seq branch: ME-T7-bulk, `sss_bulk-nXT` and standard Nextera index PCR are
  not drawn. The Methods name the bulk second-strand oligo `sss_scnXTv2`, whereas Table S1
  says `sss_bulk-nXT` is for the bulk transposon.
- Not inspected in this capped pass: the linked protocols.io Bench-top CUT&Tag
  (pre-purification steps), the review history PDF, and the Lake 2018 / Kaya-Okur 2019
  predecessors.
