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

## 6. Final library and reads

🟡 For the representative CS1 3′-end guide:

```
P5 · Nextera Read 1 · CBC16 · UMI12 · CS1 · CR1' · protospacer' · TSO' · Read 2' · i7' · P7'
```

- 🟢 Read 1 bases 1–16: cell barcode; 17–28: UMI.
- 🟡 Read 2 first traverses 30 TSO-derived bases, then bases 31–50 identify the
  protospacer, followed by the scaffold.
- 🟢 Index 1 identifies the sample in the v3 single-index workflow.

## 7. Version boundary

The chemistry shown is the v3 **CG000184** library. The v3.1 **CG000316** successor uses
dual indexing and should be modeled as a separate commercial version if its full final
construct is added; this page does not splice its index structure into the older assay.
