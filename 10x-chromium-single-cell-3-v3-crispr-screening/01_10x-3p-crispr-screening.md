# 10x Chromium 3′ gRNA capture with Feature Barcoding

This is a commercial protocol with no defining paper. Primary sources are 10x Genomics
user guide **CG000184 Rev A** and guide-RNA specification **CG000197 Rev A**. The current
v3.1 dual-index successor is **CG000316 Rev E**.

Evidence: 🟢 printed by 10x · 🟡 computed from printed sequences · 🔴 unresolved.

## 1. What it is

The 10x 3′ CRISPR Screening assay directly captures an engineered guide RNA as a Feature
Barcode. The sgRNA carries the reverse complement of Capture Sequence 1 or 2; a matching
gel-bead primer supplies the 16-base cell barcode and 12-base UMI during reverse
transcription. The guide cDNA is enriched separately from gene-expression cDNA and Read 2
reports the protospacer. 🟢

This differs from 5′ direct-capture Perturb-seq: the 10x 3′ assay changes the guide RNA
and uses a capture primer already present on the bead; Replogle 5′ capture leaves a
standard guide unchanged and spikes a guide-specific RT primer into the GEM reaction.

## 2. Sources

- 🟢 **CG000184 Rev A**: v3 reaction order, oligos, intermediate products and single-index library.
- 🟢 [CG000197 Rev A](https://cdn.10xgenomics.com/image/upload/v1660261286/support-documents/CG000197_GuideRNA_SpecificationsCompatible_withFeatureBarcodingtechnology_forCRISPRScreening_Rev-A.pdf): compatible CS1/CS2 guide designs.
- 🟢 [CG000316 Rev E](https://www.10xgenomics.com/support/universal-three-prime-gene-expression/documentation/steps/library-prep/chromium-single-cell-3-reagent-kits-user-guide-v-3-1-chemistry-dual-index-with-feature-barcoding-technology-for-crispr-screening): current v3.1 dual-index successor. This page retains the v3 CG000184 construct rather than mixing versions.
- 🟢 [10x Chromium i7 Multiplex index map](https://cdn.10xgenomics.com/raw/upload/v1655155614/support/in-line%20documents/chromium-shared-sample-indexes-plate.json): the v3 single-index kit has four 8-nt i7 oligos per well; for SI-GA-A1 these are `GGTTTACT`, `CTAAACGG`, `TCGGCGTC`, `AACCGTAA`. The [10x index-set page](https://www.10xgenomics.com/support/universal-three-prime-gene-expression/documentation/steps/library-prep/sample-index-sets-for-single-cell-3) identifies this kit for v3.

## 3. Capture sequences and guide design

🟢 Printed 5′→3′:

```
Capture Sequence 1     TTGCTAGGACCGGCCTTAAAGC
Capture Sequence 1 rc  GCTTTAAGGCCGGTCCTAGCAA
Capture Sequence 2     CCTTAGCCGCTAATAGGTGAGC
Capture Sequence 2 rc  GCTCACCTATTAGCGGCTAAGG
```

CG000197 permits either capture sequence in a scaffold hairpin or immediately before the
Pol III terminator. The schematic uses the 3′-end CS1 configuration:

```
5′- protospacer · CR1 scaffold · CS1 reverse complement · UUUUUUU -3′
```

The capture sequence is constant; the 20-base protospacer is the Feature Barcode that
identifies the perturbation. 🟢

## 4. Oligos

🟢 The CS1 gel-bead primer consists of partial Nextera Read 1N, CBC16, UMI12 and CS1.
The guide arm uses the ordinary 10x template-switch oligo after RT.

```
Feature cDNA forward  GCAGCGTCAGATGTGTATAAGAGACAG
Feature cDNA reverse  AAGCAGTGGTATCAACGCAGAG
Feature SI forward    AATGATACGGCGACCACCGAGATCTACACTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
Feature SI reverse    GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTAAGCAGTGGTATCAACGCAGAG
```

🟡 The model decomposes Feature SI forward as P5 + Nextera Read 1, and Feature SI reverse
as TruSeq Read 2 + the partial TSO primer.

## 5. Workflow

1. 🟢 Express a compatible CS1- or CS2-bearing sgRNA in the screened cells.
2. 🟢 Form GEMs. The matching capture-sequence bead primer anneals to the engineered
   guide and copies toward the protospacer.
3. 🟡 RT reaches the guide 5′ end, adds CCC and template-switches onto the 10x TSO.
4. 🟢 Feature cDNA Primers 1 enrich the short guide product between the Read 1N bead end
   and partial TSO end. Gene-expression cDNA is processed as a separate library.
5. 🟢 Feature SI PCR adds P5/full Read 1N and a TruSeq Read 2 end. Sample-index PCR adds
   P7 and the i7 index in the v3 workflow.

🟢 CG000184 Rev A detail (bounded pass 2026-10-17; txt lines): GEM-RT 53 °C 45 min,
85 °C 5 min (3282–3310); Recovery Agent + Dynabeads MyOne SILANE cleanup (2.1); cDNA
amplification with Feature cDNA Primers 1 PN-2000096, 98 °C 15 s / 63 °C 20 s / 72 °C 1 min,
11–13 total cycles (3769–3860); 0.6X SPRIselect, 150 µl supernatant kept for the CRISPR
library and pellet for 3′ GEX (3960–4000); supernatant 1.2X cleanup, 50 µl eluate
(4000–4060); step 4.1 1.0X Guide RNA cDNA Cleanup, 30 µl (5025–5060); step 4.2 Feature PCR
on 5 µl with Feature SI Primers 1 PN-2000098, 98 °C 20 s / 60 °C 5 s / 72 °C 5 s, 15 cycles
(5080–5170); 4.3 1.0X cleanup; 4.4 Sample Index PCR on 5 µl with Feature SI Primers 2
PN-2000099 + Chromium i7, 98 °C 20 s / 54 °C 30 s / 72 °C 20 s, 5 cycles (5240–5330);
4.5 0.7X/1.0X double-sided selection (5349–5380).

Published discrepancies (shown on the schematic):

- Appendix p.70 heads the Step 4.2 Feature PCR primers "Feature SI Primers 2 PN-2000099"
  while the step 4.2 reaction table uses Feature SI Primers 1 PN-2000098; the printed pair
  (P5/Read 1N + Read 2/partial TSO) is used and named PN-2000098 in the model.
- Appendix p.70 Step 4.4 product omits the 30 TSO-derived bases that its Step 4.2 product
  contains (its label row still lists "TS0"). Neither 4.4 primer removes them, so the
  model keeps them; the 4.4 drawing is read as an abbreviation.

## 6. Final library and reads

🟡 For the representative CS1 3′-end guide:

```
P5 · Nextera Read 1 · CBC16 · UMI12 · CS1 · CR1' · protospacer' · TSO' · Read 2' · i7' · P7'
```

- 🟢 Read 1 bases 1–16: cell barcode; 17–28: UMI.
- 🟡 Read 2 first traverses 30 TSO-derived bases, then bases 31–50 identify the
  protospacer, followed by the scaffold.
- 🟢 CG000184 Rev A step 5 (txt 1455–1463, 5540–5580): recommended 28/8/0/91 cycles
  (Read 1/i7/i5/Read 2); minimum Read 2 for CRISPR Screening libraries 70 bp. The model now
  acquires 91 Read 2 cycles (it previously used 50, which is below the stated minimum).
- 🟢 Index 1 identifies the sample in the v3 single-index workflow.

## 7. Version boundary

The chemistry shown is the v3 **CG000184** library. The v3.1 **CG000316** successor uses
dual indexing and should be modeled as a separate commercial version if its full final
construct is added; this page does not splice its index structure into the older assay.
