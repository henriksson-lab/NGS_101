# Modified SMART-Seq — research notes

Defining protocol: Zheng *et al.* (2022), "Single Cell Gene Expression Profiling of Embryoid Body", *Methods Mol Biol* 2520:233–259, DOI [10.1007/7651_2021_435](https://doi.org/10.1007/7651_2021_435). The author-hosted PDF (xichenlab.com) was read: Introduction, Materials 2.2 with Table 1, Methods 3.5–3.9 and Notes 1–16.

The authors describe it as a modified SMART-seq3 done in 384-well format with barcoded RT primers, pooling after cDNA amplification and 3′-end counting by UMIs. It does not replace the separate SMART-seq, SMART-seq2 or SMART-seq3 pages; its TSO is the SMART-seq2-type `AAGCAGTGGTATCAACGCAGAGTACATrGrGrG`.

## Oligos (Table 1, verbatim)

- 🟢 Barcoded RT Primer: `CTACACGACGCTCTTCCGATCT[8-bp well barcode]NNNNNNNNNNTTTTTTTTTTTTTTTTTTTTTTTTTVN` (well barcodes in Table 2).
- 🟢 TSO: `AAGCAGTGGTATCAACGCAGAGTACATrGrGrG` (rG = riboguanosine).
- 🟢 cDNA PCR Fwd / Rev.: `CTACACGACGCTCTTCCGATCT` / `AAGCAGTGGTATCAACGCAGAG`.
- 🟢 Library PCR Universal Primer: `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTC`.
- 🟢 Library PCR Index Primer: `CAAGCAGAAGACGGCATACGAGAT[8-bp plate barcode]GTCTCGTGGGCTCGG` (plate indices in Table 3; Note 11 calls it basically Illumina's Nextera XT Index PCR Primer).

🟡 Derived in the model: the universal primer is P5, a short tail (the final library shows it) and a 3′ end that overlaps the RT handle, so the final P5 end carries the full TruSeq Read 1 site. The index primer arms are the reverse complements of P7 and of the Nextera i7 (Read 2-type) Tn5 adapter end.

## Steps (Methods 3.5–3.8)

1. 🟢 Lysis plate: 2 µl lysis buffer (PEG-8000, Triton X-100, RNase inhibitor, dNTPs) plus 1 µl of 2 µM Barcoded RT Primer per well; one sorted cell per well (P24 empty control); 72 °C 10 min.
2. 🟢 RT and template switching in one program with 1 µl RT mix (Maxima H-minus, TSO, GTP, DTT): 42 °C 90 min, 10 × (50 °C 2 min, 42 °C 2 min), 85 °C 5 min.
3. 🟢 cDNA PCR in each well with KAPA HiFi HotStart ReadyMix and the cDNA PCR primers, 16 cycles.
4. 🟢 Pool all wells of one plate (plates kept separate, Note 5), purify on one QIAquick column, then **Exonuclease I** (37 °C 30 min, 80 °C 2 min) and 0.7× SPRI. Note 6: primer digestion prevents cross-well mispriming in the later PCR.
5. 🟢 Tagment 50 ng cDNA with Tn5 (TDE1, Fapon, Vazyme or in-house Tn5, Note 10) in TD buffer, 55 °C 5 min; stop with 0.2% SDS; 0.8× SPRI.
6. 🟢 Library PCR with universal and plate-index primers (KAPA HiFi; 72 °C 5 min, 98 °C 2 min, 8 cycles); 0.8× SPRI, repeated if primer leftover remains (Fig. 7).

The former model put Exonuclease I before pooling, kept the TSO-derived end inside the final 3′ library, omitted the oligo(dT) tract and the universal-primer tail bases between P5 and the RT handle. Those were repository modelling errors, now corrected.

## Read layout (3.8 step 11, Note 13)

🟢 Read 1 ≥18 cycles: 8-bp well barcode then 10-bp UMI. Index 1, 8 cycles: plate barcode. Index 2, 8 cycles, identity N/A; Note 13 calls it optional because the left side has no index. Read 2 ≥50 cycles (normally 75): cDNA.

## Open questions

- Tn5 adapter oligos are not printed; the Read 2-type end is inferred from the index primer. The fate of Tn5–Tn5 and Tn5–TSO-end fragments is not described.
- Table 2/3 barcode members were not transcribed into the model (whitelist features only).
