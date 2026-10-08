# NEBNext Ultra II DNA Library Prep — USER adaptor

## 1. What it is

This is the classic NEBNext Ultra II DNA workflow with the public non-indexed hairpin
adaptor: repair and dA-tail DNA, ligate one hairpin at each end, cut both loops open with
USER Enzyme, and complete a dual-index Illumina library by PCR.

Commercial defining sources: NEB manuals **E7103/E7645 v7.0_9/22** and **E7600
v6.0_6/24**.

## 2. Sources read

- 🟢 NEBNext Ultra II DNA Library Prep Kit for Illumina, E7103/E7645 manual v7.0_9/22.
- 🟢 NEBNext Multiplex Oligos for Illumina, Dual Index Primers Set 1, E7600 manual
  v6.0_6/24.
- 🟢 NEB USER Enzyme product documentation.
- 🟢 Illumina Adapter Sequences, document 1000000002694 v22, September 2025.

The newer indexed UMI adaptor is deliberately out of scope: it does not use USER, and NEB
now states that its correct sequence is proprietary.

## 3. End prep and ligation

🟢 End prep acts at 20 °C for 30 min and 65 °C for 30 min, repairing ends,
phosphorylating 5′ ends and adding single 3′ dA overhangs. Hairpin ligation is 20 °C for
15 min.

🟢 The E7600 adaptor is printed as
`/5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTC[dU]ACACTCTTTCCCTACACGACGCTCTTCCGATC-s-T`.
It is 65 nt: a 12-bp stem, 40-nt loop, dU at position 32 and a phosphorothioate-protected
3′ dT overhang.

🟡 One adaptor at each dA-tailed end yields a covalently closed dumbbell before USER.
`Hairpin` derives and validates the stem; `Dumbbell` derives the one closed molecular path
and the four insert–adaptor junctions.

## 4. USER opening

🟢 USER combines UDG and Endonuclease VIII. UDG removes the loop dU base and Endonuclease
VIII removes the abasic residue, leaving a one-nucleotide gap with a 3′ phosphate on the
left arm and 5′ phosphate on the right arm. Treatment is 37 °C for 15 min.

🟡 Both displayed arms come from `Hairpin.user_open_at()`; their sequences are not copied
independently from the intact adaptor.

## 5. PCR-completed library

🟢 The ligated adaptor has a deliberately truncated design, so NEB requires at least three
PCR cycles. The concrete E7600 example uses i501 (`P5 + TATAGCCT + TruSeq Read 1`) and
i701 (`P7 + CGAGTAAT + TruSeq Read 2`).

🟡 The finished top strand carries P5, i5 `TATAGCCT`, Read 1, insert, the dA junction,
the Index 1 / Read 2 side, i7 read `ATTACTCG`, and P7 reverse complement. Read 1, Index 1,
Index 2 and Read 2 primer sites are located computationally on that final molecule.
