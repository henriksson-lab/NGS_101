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
- 🟢 Adamson lab [5′ Perturb-seq Direct Capture Protocols](https://badamsonlab.squarespace.com/s/5-Perturb-seq-Direct-Capture-Protocols.pdf): exact oJR160/oJR163/oJR165 sequences, spike-in amount, PCR and size-selection details.

## 3. Oligos

🟢 Printed 5′→3′ in the laboratory protocol:

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
5. 🟢 Use 5 ng of the guide fraction in four KAPA reactions with 0.6 µM each oJR163 and
   oJR165; amplify for 12 cycles (98 °C denaturation, 70 °C anneal/extend), pool, and
   clean at 0.8×. The guide library is reported at approximately 250 bp.
6. 🟢 Spike the guide library at about 10% beside the gene-expression library and PhiX.

## 5. Final library and sequencing

🟡 The exact printed primers and the CR1 sequence produce a 251-nt model, consistent
with the protocol’s approximate 250-bp product:

```
P5 · Read 1 · CBC16 · UMI10 · switch spacer · GGG · protospacer20 · CR1 · Read 2 · i7' · P7'
```

- 🟢 Read 1 cycles 1–16: CBC; cycles 17–26: UMI.
- 🟢 Index 1: 8-cycle sample index.
- 🟢 Read 2 crosses guide-constant sequence before the protospacer; the lab protocol
  recommends **98 cycles** so the guide identity is reached. The paper’s UPR experiment
  used 26-bp Read 1, 125-bp Read 2 and 8-bp Index 1.

## 6. Boundary and open questions

- The paper also tests a cs1-modified 5′ guide and two 3′ capture-sequence designs. They
  require different molecular models and are deliberately not implied by this page.
- The method was independently paralleled by the 5′ guide-capture strategy in ECCITE-seq;
  that is a separate protocol and paper, not an alternate citation for this construct.
