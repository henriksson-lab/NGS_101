# 5′ direct-capture Perturb-seq

**Replogle JM, Norman TM, Xu A, Hussmann JA, Chen J, Cogan JZ, Meer EJ, Terry JM,
Riordan DP, Srinivas N, Fiddes IT, Arthur JG, Alvarado LJ, Pfeiffer KA, Mikkelsen TS,
Weissman JS, Adamson B.** “Combinatorial single-cell CRISPR screens by direct guide RNA
capture and targeted sequencing.” *Nature Biotechnology* 38, 954–961 (2020).
[doi:10.1038/s41587-020-0470-y](https://doi.org/10.1038/s41587-020-0470-y).

Evidence: 🟢 source text or printed oligo · 🟡 computed from those sequences · 🔴 unresolved.

## 1. What it is

The 5′ direct-capture Perturb-seq method adds one guide-specific primer to 10x Chromium
5′ reverse transcription. The primer copies a non-polyadenylated sgRNA toward its
protospacer; template switching then transfers the bead’s cell barcode and UMI onto that
guide cDNA. A separate PCR makes the short guide library, which is sequenced beside the
ordinary gene-expression library. 🟢

This page covers the paper’s **standard, unmodified sgRNA-CR1** branch using oJR160. It
does not combine the modified sgRNA-CR1cs1 branch or either 3′ direct-capture design. 🟢

## 2. Sources read

- 🟢 The defining paper and Methods in [PMC7416462](https://pmc.ncbi.nlm.nih.gov/articles/PMC7416462/): method design, 5′ capture mechanism, UPR experiment and run setup.
- 🟢 Published Supplementary Information (Springer `41587_2020_470_MOESM1_ESM.pdf`): Supplementary Note 3, "Protocol for 5' direct capture Perturb-seq of unmodified sgRNAs" (extracted lines 304–425) prints oJR160/oJR163/oJR165, the 5-pmol spike-in, the SPRI fractionation, the guide-library PCR table and the sequencing recommendation; Supplementary Notes 1–2 (lines 185–302) define CR1.
- 🟢 Published Supplementary Table 1 (PMC supplement 3 workbook): the printed guide-scaffold sequence of pBA900 (CR1 with cs2 appended) supplies the CR1 bases.
- Not read (not in the source bundle): Adamson lab [5′ Perturb-seq Direct Capture Protocols](https://badamsonlab.squarespace.com/s/5-Perturb-seq-Direct-Capture-Protocols.pdf) and the 10x Chromium Single Cell V(D)J User Guide CG000086 that Supplementary Note 3 follows. Earlier versions of this note attributed the oligo and PCR details to the lab protocol; they are now traced to Supplementary Note 3.

## 3. Oligos

🟢 Printed 5′→3′ in Supplementary Note 3 (and oJR160 also in the Methods):

```
oJR160  AAGCAGTGGTATCAACGCAGAGTACCAAGTTGATAACGGACTAGCC
oJR163  AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
oJR165  CAAGCAGAAGACGGCATACGAGATAGGAGTCCGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGAGTACCAAGTTGATAACGGACTAGCC
```

🟡 Sequence decomposition performed by the model:

- oJR160 = the 25-nt 10x RT-adapter sequence + a 21-nt guide-constant-region primer.
- oJR163 is the complete P5/TruSeq Read 1 amplification primer.
- oJR165 = P7 + the printed 8-nt example i7 (`AGGAGTCC`) + Nextera Read 2 arm + a
  nested 26-nt guide-cDNA site. Its last 21 nt equal the oJR160 guide-binding end.

## 4. Workflow

1. 🟢 Add **5 pmol oJR160** to each 68.3 µL 10x Chromium Single Cell 5′ RT master mix.
   The guide-specific primer is approximately 5% of the amount of the kit RT oligo.
2. 🟢 oJR160 anneals in the sgRNA-CR1 constant region. RT copies the constant region and
   protospacer, adds non-templated cytosines, and template-switches onto the barcoded
   bead oligo. The resulting guide cDNA carries the same CBC and UMI system as mRNA.
3. 🟢 Amplify cDNA for **11 cycles** with the 10x Non-Poly(dT) cDNA amplification primer.
4. 🟢 At the 0.6× left-side SPRI step, retain the larger gene-expression cDNA on the
   beads and take the guide-containing supernatant. Perform a 0.6×–1.2× double-sided
   selection; the guide-cDNA product is reported at about 168 bp.
5. 🟢 Use 5 ng of the guide fraction in four KAPA HiFi ReadyMix reactions with 0.6 µM
   each oJR163 and oJR165 (3 µL of 10 µM per 50 µL); 95 °C 3 min, then 12 cycles of
   98 °C 15 s / 70 °C 10 s, 72 °C 1 min; clean at 0.8× SPRI and elute in 40 µL. The
   guide library is reported at approximately 250 bp.
   **Published discrepancy:** the Methods ("Sequencing library preparation") state
   "0.6 mM oJR163 and oJR166" for the cs1 branch, whereas Supplementary Note 3 step 5
   states 0.6 uM, consistent with its reaction table. The model uses 0.6 µM; whether
   "mM" is a misprint is not stated. The schematic page shows this discrepancy.
6. 🟢 Spike the guide library at about 10% beside the gene-expression library and PhiX.

## 5. Final library and sequencing

🟡 The exact printed primers and the CR1 sequence produce a 251-nt model, consistent
with the protocol’s approximate 250-bp product:

```
P5 · Read 1 · CBC16 · UMI10 · switch spacer · GGG · protospacer20 · CR1 · Read 2 · i7' · P7'
```

- 🟢 Read 1 is 26 cycles (Methods, UPR sequencing). 🟡 The 16-nt CBC / 10-nt UMI split
  and the bead oligo's partial Read 1, switch spacer and rGrGrG follow the 10x 5′ kit
  design; they are not printed in the inspected sources and are drawn as inferred.
- 🟢 Index 1: 8-cycle sample index.
- 🟢 Read 2 crosses guide-constant sequence before the protospacer; the lab protocol
  recommends **98 cycles** so the guide identity is reached. The paper’s UPR experiment
  used 26-bp Read 1, 125-bp Read 2 and 8-bp Index 1.

🟡 CR1 bases: Supplementary Note 2 calls CR1 the parental constant region of
pU6-sgRNA EF1Alpha-puro-T2A-BFP, and Supplementary Table 1 prints pBA900's scaffold as
CR1 followed by cs2. The model takes CR1 as that printed scaffold upstream of cs2
({{= len(direct_capture.SCAFFOLD) }} nt). It has A, not the T of the Chen 2013 F+E scaffold,
at position 13 (`GAGCTAaGCTGG` as printed). An earlier model used the F+E scaffold; that
was a repository modelling error, not a publication discrepancy.

## 6. Boundary and open questions

- 🔴 The defining paper's Methods identify the 10x Poly-dT RT Primer (PN-2000007) and
  Non-Poly(dT) amplification primer (PN-220106), and the companion Supplementary Note 3
  prints the oJR oligos. The inspected sources do not print the kit bead/TSO oligo or
  Non-Poly(dT) primer sequence (the 10x CG000086 guide was not inspected; this is not a
  confirmed absence). The model carries the preamplification primers as roles and draws
  the bead oligo's partial Read 1, switch spacer and rGrGrG as INFERRED kit-design bases.
- The paper also tests a cs1-modified 5′ guide and two 3′ capture-sequence designs. They
  require different molecular models and are deliberately not implied by this page.
- The method was independently paralleled by the 5′ guide-capture strategy in ECCITE-seq;
  that is a separate protocol and paper, not an alternate citation for this construct.
