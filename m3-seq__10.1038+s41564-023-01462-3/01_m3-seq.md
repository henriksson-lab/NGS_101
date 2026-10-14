# M3-seq — research notes

Defining source: McNulty et al., 2023, [doi:10.1038/s41564-023-01462-3](https://doi.org/10.1038/s41564-023-01462-3).

## Sources read

- 🟢 PMC full text and Methods.
- 🟢 Supplementary Figure 1 workflow and Supplementary Table 3 oligo workbook.
- 🟢 The authors' `brwaang55/m3seq_scripts` repository linked by the paper, including the
  demultiplexing read structure.

## Published oligos recovered

🟢 Supplementary Table 3 prints 288 phosphorylated round-one RT primers. Each has the
structure `TruSeq Read 1 (33 nt) – N8 UMI – 13-nt BC1 – N6 random primer`; Barcode 1 uses
BC1 `AAGTGATTAGCAA`.

🟢 The same table prints the full bridge oDS025
`CGTCGTGTAGGGAAAGAGTGTGACGCTGCCGACGA/3ddC/`. Its sequence is the reverse complement
of the 10x s5 arm followed by the first 21 bases of the RT-primer handle, fixing the
splinted nick geometry rather than merely naming a ligation.

🟢 It also prints the S³ second-strand randomer oBW140, P5 primers oDS027/oBW170 and
oBW171, Read-2-only Tn5 oligos oDS029/oDS030, the T7–s7 primer oBW168, and the complete
indexed P7–s7 primer series oDS033–oDS040 and oBW176–oBW183. The page uses oBW176 as a
representative final i7 primer.

## Molecular path

🟢 Round-one random RT installs a 13-base BC1 and 8-base UMI in fixed bacteria. Round-two
10x ATAC droplet ligation installs a 16-base BC2. BC1+BC2 is the cell identity.

🟢 After second-strand synthesis the library is tagmented, PCR-tagged with a T7 promoter, transcribed, hybridized to rRNA probes and RNase-H depleted. A second RT plus PCR creates the sequenced cDNA library.

🟢 The paper specifies 26-cycle Read 1, 30-cycle i5, 8-cycle i7 and 74-cycle Read 2.
The authors' demultiplexing scripts further identify the 13-base BC1, 8-base UMI and
16-base BC2 windows.

🔴 The actual sequence assigned to any particular commercial 10x gel-bead barcode is not
printed in the paper, supplement, or linked analysis repository. The page therefore shows
its published length and molecular position as a role token, while retaining the exact
public adapter arms and bridge geometry.
