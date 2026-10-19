# Olink Explore / NGS-PEA — source notes

Defining source: Wik *et al.*, *Molecular & Cellular Proteomics* (2021),
[doi:10.1016/j.mcpro.2021.100168](https://doi.org/10.1016/j.mcpro.2021.100168).

## Evidence read

- 🟢 Full article and detailed Experimental Procedures, PMC8633680.
- 🟢 Figure 1 architecture and the PCR1, PCR2, sequencing and data-processing sections.
- 🟢 Supplementary figures and tables for controls, abundance blocks and validation.
- 🟢 Olink patent family WO2021191442A1, whose worked protocol matches the
  defining paper's abundance-block, PCR1, PCR2 and panel-pooling workflow. Its
  sequence listing publishes the P5 and P7 adapter portions of the PCR primers.
- 🟢 Two readers checked the linked Olink software landing page and current vendor
  preprocessing guide. The software page itself says its contents are password
  protected, so its archives remain inaccessible. The guide is NPX Map Preprocessing
  v6.2.2 (2026-05-19); it covers `ngs2counts` and run handling without an oligo
  table. The linked CRUK Cambridge `combine-olink` example files still need source
  review. These later materials cannot prove the 2021 probe sequences absent.

## Molecular workflow

🟢 Each protein assay uses a matched pair of antibody-conjugated oligos. One arm
contains P5 and the Read 1 sequencing-primer site; the other contains a common
reverse preamplification site. Both carry assay-specific barcodes and complementary
hybridization arms. Co-binding to the same protein raises their effective local
concentration, allowing the arms to hybridize.

🟢 The PCR1 mix (article) or extension mix (patent worked example) is added to each
incubation reaction or plate, with no wash described between the immunoreaction and
extension (JATS lines 164–168; patent text lines 4128–4135, 4994–4996). The article
says only that proximal matched pairs generate detectable amplicons (Fig. 1A legend,
line 152). 🟡 Free and non-proximal probes are therefore inferred to stay in the
reaction; the model leaves their later physical fate unresolved rather than calling
them discarded.

🟢 Samples are used undiluted and at 1:10, 1:100 and 1:1000; 0.2 µl sample is mixed
with 0.6 µl block incubation mix (80 µl incubation solution + 10 µl forward + 10 µl
reverse probes) and incubated at 4 °C for 16–24 h (PMC8633680 JATS lines 162–164).

🟢 PCR1 mix (Olink PCR1 solution, PCR1 enhancer, PCR1 Enzyme and IDT universal
primers) is added; PCR1 begins with 20 minutes at 50 °C for polymerase extension
across the paired arms, then 95 °C 5 min and 25 cycles (95 °C 30 s, 54 °C 1 min,
60 °C 1 min). Products from four concentration blocks are pooled for each sample.
PCR2 uses the same forward primer and one of 96 reverse primers to add an inline
sample index and P7 during ten cycles (95 °C 30 s, 68 °C 1 min). All PCR2 samples
of a panel are pooled into one of four libraries, and each library (50 µl) is
purified with 80 µl AMPure XP beads to remove primers, primer dimers and dNTPs,
then sized on a Bioanalyzer HS chip (JATS lines 168–172). The patent's worked
example names Pwo polymerase for extension/PCR1 (line 4135); the article names only
"PCR1 Enzyme", so the model does not assign the enzyme identity.

🟢/🟡 Apparent wording discrepancy, shown on the schematic: the article calls the
PCR1 forward primer "complementary to Illumina's P5 adapter sequence" (JATS line
168), whereas the patent says the forward common primer "comprises" the P5 adapter
(SEQ ID NO: 1; lines 4135–4171). The model uses the P5-sense identity, which the
article's own probe design (forward probe carries P5 and Rd1SP) requires.

🟢 The finished amplicon is 148 bp. Its published element order is P5, Rd1SP,
forward assay barcode, hybridization-derived junction, reverse assay barcode,
sample index and P7.

🟢 WO2021191442A1 SEQ ID 1 is the 20-nt P5 adapter portion
`AATGATACGGCGACCACCGA`; the patent states that the PCR1 forward common primer
contains it and the same primer is reused in PCR2. SEQ ID 2 is the 24-nt P7
adapter portion `CAAGCAGAAGACGGCATACGAGAT`; each PCR2 reverse primer contains
it plus a sample-index barcode. These are published portions, not the complete
primer sequences.

🟡 The rest of the exact standard Illumina P5/Rd1SP end shown on the schematic
comes from Illumina's adapter documentation. The combined end is 58 nt. Together
with the patent's 24-nt P7 adapter, the published 148-bp total leaves exactly 66 nt
between Rd1SP and P7. This agrees with the published 66-cycle read covering the inline
assay and sample identifiers; it is a consistency check, not proof of the bases, and
the article does not state where Read 1 starts or ends. The individual lengths
within those 66 nt have not yet been verified from an original oligo table.

🔴 The original article, its four supplements and the checked patent do not supply
the complete probe oligos, universal PCR primers, assay-barcode whitelist,
hybridization arms or sample-index sequences in the inspected sections. This is a
provisional source gap: the software archive is password protected, the checked
preprocessing guide is from 2026, and third-party example counts have not yet
received an independent second source check.
The schematic therefore shows these regions by identity and order without inventing
bases or individual lengths. No raw-read accession was found in the article's data
availability statement; this absence claim also remains under review.

## Additional public material checked

- [Olink software downloads](https://olink.com/software/download)
- [WO2021191442A1 — Method for detecting analytes of varying abundance](https://patents.google.com/patent/WO2021191442A1/en)
- [Olink Explore preprocessing technical information](https://7074596.fs1.hubspotusercontent-na1.net/hubfs/7074596/01-User%20Manuals%20for%20website/1344-Olink-explore-preprocessing.pdf)
- [CRUK Cambridge `combine-olink`](https://github.com/crukci-bioinformatics/combine-olink)

## Sequencing

🟢 Each of four panel libraries is loaded in one NovaSeq 6000 XP lane. A single
66-cycle Read 1 covers both assay barcodes and the inline sample index. There are no
separate index reads: exact barcode combinations are counted to estimate protein
abundance after control normalization.

## Review trail

- 2026-10-17 bounded pass (round 2026-10-17-breadth-first-14): re-read JATS lines
  150–176 and patent lines 1746–1753, 4135–4195; added the 4 °C incubation
  condition, PCR1-mix extension catalyst, panel pooling and AMPure XP cleanup with
  its discarded supernatant, and the forward-primer wording discrepancy. Not
  re-read: supplements mmc1–mmc4, the 2026 preprocessing guide, `combine-olink`.
- 2026-10-19 bounded migration check (round 2026-10-19-migration-final): removed the
  unused shared `batch_ngs` imports (the graph was already native); re-read JATS
  lines 152, 158–176 and patent lines 1753, 4128–4135, 4996. Unbound probes changed
  from a physical discard to an unresolved terminal (no wash step); the linker and
  Rd1SP bases between patent SEQ ID 1 and the inline region are now marked inferred
  on the page, matching the 🟡 label above.
- 2026-10-19 independent migration review (batch14_independent): re-read JATS lines
  152, 164–176 and patent lines 1753, 4128–4135, 4946–4996. The no-wash statement is
  now 🟢 only for the described mix addition, with probe carry-over 🟡; the 58 + 66
  + 24 = 148 arithmetic was confirmed as a consistency check, and the page caption no
  longer calls the derived 66-nt inline length published.
