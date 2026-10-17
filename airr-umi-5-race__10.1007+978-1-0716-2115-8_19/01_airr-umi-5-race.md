# Bulk UMI 5′-RACE AIRR-seq

## Source read

🟢 Gupta et al., *Bulk Sequencing from mRNA with UMI for Evaluation of B-Cell Isotype and Clonal Evolution: A Method by the AIRR Community*, in Langerak (ed.), *Immunogenetics: Methods and Protocols*, Methods in Molecular Biology vol. 2453 (2022), chapter 19, [doi:10.1007/978-1-0716-2115-8_19](https://doi.org/10.1007/978-1-0716-2115-8_19). Open access (CC BY) as NCBI Bookshelf NBK586949; the local copy is `ch19.pdf`/`ch19.nxml` from the Bookshelf open-access package of the whole book (NBK586966), identity checked by the printed DOI and title.

Sections read (2026-10-17 bounded pass, `NBK586949_ch19.pdf.txt`): Introduction and Fig. 1–2 captions; 2.1 reagents; 2.2 index tables; 3.3 first-strand synthesis; Table 1; 3.4–3.5 PCR1/PCR2; 3.6 bead purification; 3.7–3.8 validation, pooling and sequencing; Fig. 5 (rendered page 359); 3.10 steps 2 and 5–6 (read layout, UMI extraction); Notes 1–14. The kit is the Takara SMARTer Human BCR IgG IgM H/K/L Profiling Kit. The Takara kit user manual was not read.

## Molecular path

🟢 dT Primer and BCR enhancer are heated with total RNA at 72 °C for 3 min, chilled, and oligo(dT)-primed SMARTScribe RT runs at 42 °C for 90 min and then 70 °C for 10 s. The RT adds non-templated nucleotides; the SMART UMI Oligo anneals to them, templates a PCR handle, and gives each cDNA a 12-nt UMI.

🟢 PCR1 (PrimeSTAR GXL, 18 cycles) pairs hBCR PCR1 Universal Forward, which anneals through the SMART UMI Oligo sequence, with one IgG, IgM, IgK or IgL constant-region reverse primer per reaction. PCR2 is semi-nested (16–25 cycles, Table 1): hBCR PCR2 Universal Forward 1–12 add P7/i7, and hBCR PCR2 Ig reverse 1–4 anneal in the constant region and add P5/i5.

🟢 Index sequences are printed (sections 2.2.2–2.2.3): i7 D701–D712 on the forward primers and i5 D501–D504 on the reverse primers, the i5 values given “as read on a MiSeq instrument”. The model holds both lists as source transcriptions and keeps the library index windows as placeholders.

🟢 Purification uses NucleoMag beads: 25 µL of beads, keep the supernatant and discard the beads; 10 µL more beads, discard the supernatant; two 80% ethanol washes; elute in 17 µL. Expected libraries run as a broad 500–1200 bp peak.

## Read layout and uncertainty

🟢 MiSeq Reagent Kit v3, 2 × 300, with 10% PhiX. Read 1 starts in the constant region and Read 2 starts upstream of the V region. Read 2 begins with a 12-nt UMI followed by a 7-base linker plus template switch (`MaskPrimers.py extract --start 12 --len 7`), and the first 19 nt may be trimmed when the UMI is not used.

🟡 Fig. 5 labels the linker as 4 bp, followed by an unsized “XXXXX” block, then the leader. Taken with the text's 7-base total, this suggests a 4-base linker plus a 3-base template-switch block. The example read annotation `PRIMER=GTACGGG` fits that split, but it is one example read, not a printed oligo. The model keeps only the published 7-base total. Fig. 5 also draws each adapter end as about 70 bp and the C-region primer plus 7N as 27–30 bp.

🔴 The chapter says “Primers from the vendor are not available”. The dT primer, SMART UMI Oligo, PCR1/PCR2 primer bodies and adapter arms are not printed in the chapter, so the model uses named role placeholders and claims no primer-binding geometry. This is a bounded single-reader finding and not an accepted absence claim. The Takara kit user manual and the Immcantation constant-region primer files are named, unread follow-ups.
