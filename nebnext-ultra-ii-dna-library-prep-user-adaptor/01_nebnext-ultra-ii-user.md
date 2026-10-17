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

The indexed UMI adaptor is a separate branch: the main manual specifies that it does
not require USER. This page selects the non-indexed E7600 hairpin and its Set 1 PCR primers;
it does not assign that hairpin structure to arbitrary third-party or indexed adaptors.

## 3. End prep and ligation

🟢 End Prep (Ultra II End Prep Reaction Buffer + Enzyme Mix) runs 20 °C for 30 min and
65 °C for 30 min (E7103/E7645 Step 1.3). The E7103/E7645 manual does not describe the
mix's activities; E7600 p. 2 says the adaptor ligates to "end-repaired, dA-tailed DNA".
🟡 5′ phosphorylation of the insert is inferred, not printed in the manuals read.
🟢 Ligation adds the adaptor, Ultra II Ligation Master Mix and Ligation Enhancer directly to
the End Prep mixture, 20 °C for 15 min (Steps 2.2–2.4).

🟢 The E7600 adaptor is printed as
`/5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTC[dU]ACACTCTTTCCCTACACGACGCTCTTCCGATC-s-T`.
It is 65 nt: a 12-bp stem, 40-nt loop, dU at position 32 and a phosphorothioate-protected
3′ dT overhang.

🟡 One adaptor at each dA-tailed end yields a covalently closed dumbbell before USER.
`Hairpin` derives and validates the stem; `Dumbbell` derives the one closed molecular path
and the four insert–adaptor junctions.

## 4. USER opening

🟢 USER Enzyme (3 µl) is added to the ligation mixture and incubated at 37 °C for 15 min
(E7103/E7645 Steps 2.5–2.6); these steps are only for the non-indexed NEBNext Adaptor, and
the indexed UMI adaptor needs no USER. E7600 p. 2 states that USER is "a combination of
UDG and Endo VIII" that removes the loop U "to open up the loop and make it available as
a substrate for PCR". The gap-end chemistry (3′ phosphate on the left arm, 5′ phosphate on
the right arm, one-nucleotide gap) comes from the USER Enzyme product documentation, which
was not re-read in the 2026-10-17 pass.

🟡 Both displayed arms come from `Hairpin.user_open_at()`; their sequences are not copied
independently from the intact adaptor.

## 5. PCR-completed library

🟢 The ligated adaptor has a deliberately truncated design, so NEB requires at least three
PCR cycles. The concrete E7600 example uses i501 (`P5 + TATAGCCT + TruSeq Read 1`) and
i701 (`P7 + CGAGTAAT + TruSeq Read 2`).

🟡 The finished top strand carries P5, i5 `TATAGCCT`, Read 1, insert, the dA junction,
the Index 1 / Read 2 side, i7 read `ATTACTCG`, and P7 reverse complement. Read 1, Index 1,
Index 2 and Read 2 primer sites are located computationally on that final molecule.


## 6. Linked FAQ and source correction (16 October 2026)

The [NEB third-party adaptor FAQ](https://www.neb.com/faqs/can-i-use-this-nebnext-kit-with-adaptors-and-primers-from-other-vendors-than-neb)
says other-vendor T-overhang adaptors may be compatible, with optimization recommended.
Its E7395 DNA and E7416 RNA links concern separate full-length indexed/UMI-adaptor
branches; compatibility does not make those adaptors USER-cleavable hairpins.

The original E7600 v6.0_6/24 **p. 3** and E7780 **p. 3** both print the full modified
NEBNext hairpin oligo. A former repository source-check exception incorrectly claimed
these manuals omitted the sequence; that exception has been removed. This was a
repository error, not a vendor publication error or demonstrated sequence withdrawal.
The adaptor and illustrative i501/i701 primer states are now checked against the actual
E7600 original, including phosphate, dU and terminal phosphorothioate notation.

## 7. Cleanups and PCR (bounded pass, 17 October 2026)

🟢 E7103/E7645 v7.0_9/22: after USER, inputs > 50 ng are size-selected with about 0.4X
beads (keep supernatant) then 0.2X (keep beads) (Section 3A); inputs ≤ 50 ng get a 0.9X
cleanup (Section 3B); both elute 17 µl and take 15 µl into PCR. PCR uses Ultra II Q5
Master Mix with one i7/Index primer and one i5 (or Universal, for single-index kits) primer:
98 °C 30 s; 3–15 cycles of 98 °C 10 s and 65 °C 75 s; 65 °C 5 min (Step 4.1.3, Tables
4.1/4.2). A 0.9X cleanup follows (Section 5). The graph draws the ≤ 50 ng 0.9X branch and
names the size-selection alternative in a page note.

Open: the USER gap-end chemistry and Illumina adapter document were not re-read; the
End Prep mix composition is not printed.

