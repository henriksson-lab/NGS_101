# Illumina 16S metagenomic sequencing

## 1. What it is

Illumina's demonstrated 16S workflow is a reusable two-stage amplicon design: PCR1 amplifies the bacterial V3–V4 locus with locus-specific primers carrying Nextera overhangs, and PCR2 uses those overhangs to add P5/P7, i5/i7 indexes and the standard sequencing-primer sites.

## 2. Authoritative source

🟢 Illumina, *16S Metagenomic Sequencing Library Preparation*, document 15044223 B. The guide gives both full PCR1 oligos, identifies the Klindworth V3–V4 primers, specifies a limited-cycle index PCR using Nextera XT dual indexes, and diagrams the final paired-end library.

🟢 The same guide's “Dual Indexing Principle” table (p. 23) gives all 12
N701–N712 i7 sequences and all 8 S501–S508 i5 sequences for the named
FC-131-1002 kit; its FC-131-1001 subset uses N701–N706 and S501–S504.
Thus the index families are published in the defining source. The guide
does not choose one pair for a particular sample, so the library schematic
retains variable index segments.

🟢 The linked Klindworth et al. 2013 original (doi:10.1093/nar/gks808,
PMC3592464 JATS line 84) independently prints the selected bacterial
S-D-Bact-0341-b-S-17 `CCTACGGGNGGCWGCAG` and
S-D-Bact-0785-a-A-21 `GACTACHVGGGTATCTAATCC` primer pair. These are exactly
the 3′ locus-specific portions of the two Illumina PCR1 oligos below. The source
does not assign a Nextera XT index pair to any Illumina-guide sample.

## 3. Oligos

🟢 Forward PCR1 primer, 5′→3′:

`TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG`

🟢 Reverse PCR1 primer, 5′→3′:

`GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC`

The 5′ portions are universal Nextera overhangs; the 3′ portions bind the 16S locus. Degenerate IUPAC bases are copied verbatim from the guide.

## 4. Construct logic

PCR1 copies the V3–V4 interval and places a different universal overhang on each end. After bead cleanup, PCR2 primers anneal to those overhangs and add the two flow-cell arms plus combinatorial i5/i7 sample indexes. The source-known overhangs supply the Nextera-compatible Read1/Read2 and Index1 sites. MiSeq Index2 uses the effective forward-strand P5 boundary described below.

## 5. Scope

🟢 Illumina explicitly states that the same tailed-primer design can target other loci. This page nevertheless remains the named V3–V4 protocol; the fungal ITS page documents a materially different EMP one-stage construct and custom primers.

## 6. Source-led model correction (2026-10-16)

The former generic two-stage endpoint discarded both gene-specific primer regions
from the declared PCR1 product. That was a repository modeling error, not a
publication error. The local endpoint now extends **the actual PCR1 product**:

```
5′-inferred P5 · i5 · forward overhang · forward locus primer · [variable biological interval] · reverse-locus reverse complement · reverse-overhang reverse complement · i7_RC · inferred P7_RC-3′
```

Both complete PCR1 primers are source-positive (guide pp. 3–4, extracted
lines110–155). Their IUPAC symbols describe nucleotide classes, not an unknown
primer length. The model computes the full IUPAC complement, including H/V↔D/B,
and declares their actual lengths. Read1 begins with the 17-base forward locus
region; Read2 begins with the 21-base reverse locus region. The biological interval
has no selected sequence or fixed length. Consequently neither read is predicted
to reach the opposite primer at a coordinate fabricated from a short placeholder.

The source guide names KAPA HiFi HotStart ReadyMix in both PCRs, 25 amplicon-PCR
cycles and 8 index-PCR cycles. It specifies AMPure cleanup after each PCR. After
PCR1 cleanup, 5µl of the 50µl eluate enters index PCR; the remaining45µl is not used
in this protocol and may be stored. The graph preserves that unused aliquot as a
separate possible stored outcome, rather than claiming it was destroyed or
molecularly rejected (guide pp.9–13, lines380–557).

PCR1 binding geometry is a model gap: no native community genomic template or
its binding coordinates has been selected. Checking the tailed PCR1 product as
though it were that native input would not prove PCR1. The full source primer
strands are still explicit graph inputs.

**INFERRED — PCR2 full strands and caps.** The local constructor defaults to
inferred canonical P5/i5/S5 and P7/i7/S7 primers; no particular index member is
selected. Its terminal14/15-base S5/S7 arms are checked against the actual PCR1
overhangs. This verifies the declared conditional geometry, not the identity of an
unassigned historical full-primer original. The defining guide positively prints
the index-family values; that does not independently establish every complete
ordered primer strand for a selected kit generation.

The source's paired300-base MiSeq v3 run is stated in the guide overview
(lines77–79), with eight-base dual indexes (lines1238–1294). Index2 follows the
MiSeq forward-strand mechanism. The independently read manufacturer *Indexed
Sequencing Overview*15057455v08 (November2020), p.4, describes graftedP5 extension
through seven non-imaged chemistry cycles before the imaged index read. The model
uses canonicalP5 as an **effective post-extension read-start anchor**, not as an
exact physical29-mer grafted primer. The seven dark cycles are additional
chemistry, not extra reported i5 bases. This contemporary comparator does not
assign the authors' historical instrument-recipe revision.

## 7. Remaining source identity and scope

A prior checklist row described a desired “historical Nextera XT full-primer
ordering document” without identifying an actual original document. It is retained
as a historical identity question, rather than being counted as an inaccessible
published disclosure. Full-primer generation and a chosen experimental sample
remain model assignments, not evidence that the authors omitted information.

The actual Klindworth article and linked-source inventory identify supplementary
Tables archives `supp_gks808_Tables.zip` and `supp_gks808_gks808-tables-xlsx.zip`,
plus the original primer-evaluation archive at `arb-silva.de/download/archive/primer_evaluation`.
These are concrete outstanding linked-source obligations; their contents have not
been excluded from relevance based on titles. Figure, script and FASTA archive
scope also remains bounded by the actual inspections recorded in the source review.
No large reference sequences or count datasets were fetched in this model-fix batch.

The review remains in progress. No missing-information claim is accepted from the
unselected sample, inferred caps, unassembled genomic template, or unchecked
linked originals.
