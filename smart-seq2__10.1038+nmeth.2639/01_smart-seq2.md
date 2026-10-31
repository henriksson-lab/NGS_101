# SMART-seq2 — oligo source review

Picelli S, Björklund ÅK, Faridani OR, Sagasser S, Winberg G, Sandberg R. "Smart-seq2 for
sensitive full-length transcriptome profiling in single cells." *Nat Methods*
2013;**10**(11):1096–1098. doi:[10.1038/nmeth.2639](https://doi.org/10.1038/nmeth.2639)

> **Evidence marking.** 🟢 verbatim from a source read · 🟡 derived or inferred ·
> 🔴 not published in any source read here.

This note exists to settle one question: where does every base on the SMART-seq2 schematic
come from? The short answer is that all of them are transcriptions of primary documents,
with one two-base exception that the page now draws as inferred.

## Sources read

- 🟢 **Picelli et al. 2013 article, publisher PDF** (`nmeth.2639.pdf`, printed
  doi:10.1038/nmeth.2639, 5 pp. incl. Online Methods; hand-downloaded by the maintainer,
  identity checked against printed DOI/title/authors 2026-10-17). Read in the `.txt`:
  main text (lines 187–330), Online Methods "Single-cell cDNA isolation" through
  "Tagmentation reaction and final PCR amplification" (lines 600–700) and "Read
  alignments" (lines 738–745).
- 🟢 **Picelli et al. 2014 *Nature Protocols*** (`nprot.2014.006.pdf`, printed
  doi:10.1038/nprot.2014.006; hand-downloaded): Materials list (lines 288–296) and Reagent
  Setup oligos (lines 325–348) only. The Procedure (Nextera XT back end) was not read; it
  is a later version and does not override the 2013 values.
- 🟢 **Picelli et al. 2013 Supplementary Table 2** (Springer
  `41592_2013_BFnmeth2639_MOESM301_ESM.xlsx`), "List of all template switching
  oligonucleotides tested" — thirteen TSOs with their 5' and 3' modifications.
- 🟢 **Picelli et al. 2013 Supplementary Table 4** (`…MOESM303…`), "Detailing variants of
  the Smart-seq2 protocol" — the row named `Smart-seq2` uses TSO **`rGrG+G`**, 1 µl of
  10 µM, no bead purification before pre-amplification, KAPA HiFi HotStart ReadyMix.
- 🟢 **Picelli et al. 2013 Supplementary Tables 1 and 3** (`…MOESM300…`, `…MOESM302…`) —
  per-library reaction tables. They identify the oligo-dT only by its commercial name,
  *SMARTer dT30VN*, and the PCR primer only as the kit's *ISPCR primers*. The single-cell
  rows using `rGrG+G` with KAPA HiFi HS (Supplementary Table 3, worksheet A) list RT enzyme
  `SSRTII`, MgCl2 column `12` mM (SMARTer kit rows: 6), betaine `1` M and no purification
  after RT.
- 🟢 **Picelli et al. 2013 Supplementary Table 5** (`…MOESM304…`), cost of "this paper
  (SMART-SEQ2+NEXTERA)": SuperScript II (Invitrogen 18064-014), betaine, KAPA HiFi HotStart
  ReadyMix, Zymo DNA Clean & Concentrator-5, **Nextera DNA sample preparation kit
  FC-121-1031**, **Nextera 24-index kit FC-121-1011**, **TruSeq Dual Index sequencing primer
  kit for Single Read runs FC-121-1003**, AMPure XP, MgCl2, SMARTER dT30VN and ISPCR primers
  (Biomers), rGrG+G TSO (Exiqon).
- 🟢 **Hagemann-Jensen et al. 2020 Supplementary Table 1**, condition 1 ("Smartseq2 —
  Smart-seq2 as of Nature Protocol") — prints that version's oligo-dT
  `AAGCAGTGGTATCAACGCAGAGTAC(T)30VN`, TSO `AAGCAGTGGTATCAACGCAGAGTACATrGrG+G` and
  forward PCR oligo `AAGCAGTGGTATCAACGCAGAGT` in full.
- 🟢 **Hahaut et al. 2022 Supplementary Note 2** — its *Smart-seq2* block prints
  `Smart dT30VN` = `/5Biosg/AAGCAGTGGTATCAACGCAGAGTAC(T)30VN` and `ISPCR` =
  `/5Biosg/AAGCAGTGGTATCAACGCAGAGT`, both credited to Picelli et al. 2013.
- 🟢 **Ramsköld et al. 2012** (PMC3467340) — the identical SMART-seq CDS primer
  `AAGCAGTGGTATCAACGCAGAGTACT(30)VN`, which is where this oligo-dT comes from.
- 🟢 **Illumina customer sequence letter, printed 24 May 2012**, public mirrored
  vendor original (`source_pdf_fbf4f1d1ac38acc3ee4b.pdf`): printed pages1–3,
  text10 and86–183. Printed page3 names the same FC-121-1031/1030 and
  FC-121-1011/1012 product families, prints both transferred adapters and both
  full index-PCR primer orders, and lists eight-base N5xx/N7xx index members.
  The rendered original confirms the adapter-versus-sample-sheet columns.
- 🟢 **Illumina Using Nextera DNA Sample Preparation Kits and Enabling Dual
  Indexing, part15029964 RevC**, printed©2012
  (`illumina_15029964_C_Nextera_DNA_Sample_Preparation_Insert.pdf`): both pages,
  text34–68 and123–136. Its catalogue IDs match Supplementary Table5, and it
  specifies eight index cycles on each side for the named kit's standard workflow.
  The rendered second page was inspected too.
- 🟢 The above2012 sequence letter directly supports ME, s5, s7, P5, P7 and
  the Nextera N5xx/N7xx primer component orders already imported through
  `lib/nextera.py` and `lib/illumina.py`; a modern sequence listing is no longer
  the sole evidence for these source-era components.

The full source-era guide15027987 RevB and the complete PPC/NPM reagent compositions
remain unassembled. The2012 vendor originals now independently establish the named
kit's adapter/primer architecture and eight-base index widths; those positive findings
supersede the earlier index-width deferral. The actual index pair assigned to each
study well remains an open study-specific assignment, not a claim of absent disclosure.
The2013 Online Methods print the oligo-dT30VN and ISPCR bases and say both were
ordered from Biomers.net, so the earlier question about the commercial *SMARTer dT30VN*
identity no longer affects the drawn oligos.

This review supersedes an earlier provenance note that credited part of the exact oligo
bases to a secondary drawing of the SMART-seq family. That is not where they come from,
and that drawing is not evidence for anything on these pages.

## Oligos

```
oligo-dTVN   5'- AAGCAGTGGTATCAACGCAGAGT AC TTTTT...TTTTT VN -3'   (30 T)   🟢 2013 Online Methods
TSO          5'- AAGCAGTGGTATCAACGCAGAGT ACAT rGrG+G -3'                   🟢 2014 / 🟡 vs 2013 (see below)
ISPCR        5'- AAGCAGTGGTATCAACGCAGAGT -3'                               🟢 2013 Online Methods
```

`AAGCAGTGGTATCAACGCAGAGT` is the 23-nt **SMART / ISPCR handle**, shared by the oligo-dT,
the TSO and the PCR primer — which is the whole architectural point of the method. `rG`
is a ribonucleotide and `+G` the locked nucleic acid that is SMART-seq2's headline
improvement over SMART-seq.

### The one disagreement between sources: AC or ACAT?

🟡 Two published forms of the TSO linker exist, and they differ by two bases:

| Form | Printed in | TSO |
|---|---|---|
| `AC` | Picelli et al. 2013 Supplementary Table 2, row `rGrG+G` — the row Supplementary Table 4 names as *Smart-seq2*'s TSO | `AAGCAGTGGTATCAACGCAGAGT` **`AC`** `rGrG+G` |
| `ACAT` | Picelli et al. 2014 *Nature Protocols* Reagent Setup (line 331); also Hagemann-Jensen et al. 2020 Supplementary Table 1, condition 1 | `AAGCAGTGGTATCAACGCAGAGT` **`ACAT`** `rGrG+G` |

The 2013 Online Methods do not print the TSO: they give "1 µl TSO (10 µM, the complete
list of the oligos can be found in Supplementary Table 1)" (lines 639–641), but Supplementary
Table 1 lists libraries by TSO name and the sequences are in Supplementary Table 2. That
is a table-pointer slip in the article, recorded on the page with the TSO discrepancy.

Supplementary Table 2 is internally consistent about this: three of its thirteen oligos
(`C6 Amino`, `rGrGrG`, `rG3p`) carry `ACAT`, and the other ten, `rGrG+G` among them, carry
`AC`. So the defining paper's own chosen oligo has the shorter linker, and the longer one
belongs to the *Nature Protocols* version of the protocol published a year later.

The model keeps the longer, `ACAT` form, because that is what the protocol as practised
uses and what `smart-seq2__10.1038+nmeth.2639/tools/smartseq2.TSO_LINKER` supplies to the SMART-seq2 derivatives
(scM&T-seq, scMT-seq, scNMT-seq, scCAT-seq). The page therefore draws **`AC` as
established and the trailing `AT` as inferred**, with all three required signals: a
preamble caveat naming both forms, an `INFERRED — …` caption on the template-switch step,
and `<inf>` around the two bases in the oligo listing and the step panels. No base is
presented as published by the defining paper when it is not.

🟡 This two-base question has no consequence for the final library: the handle region the
TSO contributes is internal to the amplified cDNA and is tagmented away, so it does not
appear in the sequenced construct at all. It matters only for ordering the oligo.

## What the page asserts, by role

| Role | Status | Where it comes from |
|---|---|---|
| oligo-dT (dT30VN with SMART handle) | 🟢 | Picelli 2013 Online Methods (line 632); Picelli 2014 Reagent Setup; Ramsköld 2012 CDS primer |
| TSO handle + `AC` + `rGrG+G` | 🟢 | Picelli 2013 Suppl. Tables 2 and 4 |
| TSO linker's trailing `AT` | 🟡 | Picelli 2014 Nature Protocols Reagent Setup (later version); drawn inferred for the 2013 method |
| ISPCR primer | 🟢 | Picelli 2013 Online Methods (line 664); Picelli 2014 Reagent Setup |
| Nextera ME / s5 / s7, read and index primer entry | 🟢 | Illumina24May2012 vendor sequence letter, printedp3, exact source-kit families; imported through `lib/nextera.py` |
| P5 / P7 flowcell ends | 🟢 | Illumina24May2012 vendor sequence letter, printedp3 index-PCR orders; imported through `lib/illumina.py` |
| i5 / i7 index widths and members | 🟢 width; 🟡 sample assignment open | Eight-base N5xx/N7xx members and exact kit IDs verified in the2012 vendor originals. Generic eight-base windows do not assign the actual study-well pair. |
| Read plan | 🟢 / 🟡 | Online Methods: "sequenced with Nextera dual indexes (i7+i5) on Illumina HiSeq 2000, giving 43bp reads after demultiplexing" 🟢; no Read 2 🟡 from the single-read primer kit (Supplementary Table 5) and the all-SINGLE deposit (ENA PRJNA213629); index cycle counts not printed |

## Open roles

- 🔴 The chosen i5/i7 index pair per well is not established and is not drawn.
- 🟡 The TSO linker's trailing `AT` is printed by the 2014 *Nature Protocols* version
  (now read) but not by the 2013 paper, whose selected `rGrG+G` row carries `AC`.
- 🟡 Complete PPC/NPM compositions and the full source-era preparation guide remain
  unassembled. The2012 source-kit adapter/index orders and index widths are now verified.

## Chemistry and steps

The chemistry and final library follow SMART-seq's, with the LNA TSO, betaine and higher
MgCl₂ as the optimisations. 🟢 From the 2013 Online Methods (lines 600–700):

1. Lysis in 2 µl 0.2% Triton X-100 + 2 U/µl RNase inhibitor; add 1 µl 10 µM oligo-dT30VN
   and 1 µl 10 mM dNTPs; 72 °C 3 min, ice.
2. RT (10 µl): 0.5 µl SuperScript II, 0.25 µl RNase inhibitor, 2 µl 5× first-strand
   buffer, 0.25 µl 100 mM DTT, 2 µl 5 M betaine, 0.9 µl 100 mM MgCl₂, 1 µl 10 µM TSO;
   42 °C 90 min, 10 × (50 °C 2 min, 42 °C 2 min), 70 °C 15 min.
3. Pre-amplification without purification: 25 µl KAPA HiFi HotStart ReadyMix, 1 µl 10 µM
   ISPCR, 14 µl water; 98 °C 3 min, 18 cycles for cells (98 °C 15 s, 67 °C 20 s, 72 °C
   6 min), 72 °C 5 min. AMPure XP 1:1, 15 µl EB. Expected 1.5–2.0 kb average.
4. Tagmentation: 5 ng cDNA, 25 µl 2× Tagment DNA Buffer + 5 µl Tagment DNA Enzyme in
   50 µl, 55 °C 5 min; Zymo DNA Clean & Concentrator-5, 20 µl RSB.
5. Enrichment PCR on the whole eluate with 15 µl NPM, 5 µl N7xx, 5 µl N5xx, 5 µl PPC:
   72 °C 3 min (the gap fill), 98 °C 30 s, 5 × (98 °C 10 s, 63 °C 30 s, 72 °C 3 min).
   AMPure XP 1:1; Qubit; 2 nM pools, 10 pmol on HiSeq 2000.

🟡 The 72 °C 3-min opening of the enrichment PCR is drawn as the Nextera gap-fill
extension; the PPC outer primers are a declared input without printed bases.

## Bounded pass 2026-10-17 (round breadth-first-18)

Read: Supplementary Tables 2 (lines 3–16), 3 (worksheet A header and the reaction
columns of all rows), 4 (lines 3–9) and 5 (lines 2–35); Supplementary Figures PDF item
list (lines 1–40). Scraper hits checked: Supplementary Table 2 rows; Table 5 TSO line.
The article, its Online Methods and the Nature Protocols companion remain unread.

## Bounded pass 2026-10-17 (round deferred-followup-02)

Recovery routes for the article and the *Nature Protocols* companion: OpenAlex
(`oa_status` closed, no repository full text), Unpaywall (no OA location, no repository
copy) and Europe PMC (PMID 24056875 / 24385147, no PMCID) all return closed; the earlier
Sandberg-lab and publisher-PDF attempts are recorded in the review file. Both remain
manual (https://doi.org/10.1038/nmeth.2639, https://doi.org/10.1038/nprot.2014.006).

🟢 GEO GSE49321 (all 76 GSM records): library construction is stated only as "according
to the optimized Smart-Seq2 protocol, variants of Smart-Seq2 and ... standard Smart-Seq"
with SMARTer kit lysis buffer; 44 samples are labelled `library prep protocol: SMART-Seq2`,
the others variants 1–4 or SMARTer; instrument Illumina HiSeq 2000. 🟢 ENA read-run report
for PRJNA213629 (SRP028301): all 76 runs are `SINGLE` layout, HiSeq 2000, base count /
read count = 43. 🟡 This supports the single-read plan already drawn; the 43-base
deposited length is not drawn as the acquired cycle count because trimming before
deposit is not excluded by these records.

## Bounded pass 2026-10-17 (round deferred-followup-03)

Read the hand-downloaded defining article (Online Methods) and the *Nature Protocols*
Reagent Setup. Graph now names the RT program, KAPA pre-amplification program, 1:1 AMPure
cleanup, tagmentation conditions, Zymo DCC-5 cleanup, NPM/PPC enrichment PCR (5 cycles)
and final AMPure cleanup; ISPCR has an exact state; 43-base processed/deposited reads are reported; acquired Read 1 cycles remain unspecified (corrected during critical re-review). The TSO
discrepancy now cites the 2014 Reagent Setup directly and notes the article's
Supplementary Table 1 pointer slip.


## Critical re-review, 2026-10-17

43-base demultiplexed/deposited reads do not establish acquired cycles. Removed the acquired R1=43 assignment. Same-end tagmentation classes have unresolved amplification/recovery rather than physical discard; the end-class validator does not establish those fates.

Inspected scope: nmeth.2639.pdf.txt Online Methods, Read alignments (740–748); existing note ENA run report and processed-read distinction. This bounded correction does not complete the linked-source or two-reader absence audit.


## Source-era Nextera gap closure, 31 October2026

🟢 Fetched and validated two vendor originals through registered extra sources into
this protocol's external source bundle. The bundle scraper was rerun, including an
exact search for the transferred s5+ME adapter. The source-line pointer was checked
in the2012 letter; upstream schematic hits were excluded. The letter's printed date,
full component orders and index-table column relationships were also checked against
rendered original pages. A source-guide filename or PDF metadata date was not treated
as printed identity.

The evidence chain is specific: the2013 Supplementary Table5 rows27–29 identify
FC-121-1031, FC-121-1011 and FC-121-1003; the2012 RevC insert lists those same IDs;
the2012 sequence letter's printedp3 explicitly assigns the transferred adapter and
index-PCR orders to those kit families. This closes the previous source-era adapter
identity and index-width gap without substituting a2025 catalogue or NexteraXT recipe.
The modern canonical sequences used in the model match these contemporaneous orders.

🟢 The original letter's N501 row prints adapteri5 `TAGATCGC` and the same string for
its source-era sample sheet. N701 prints adapteri7 `TCGCCTTA` and sample-sheeti7
`TAAGGCGA`; the latter is its reverse complement. A named N501/N701 pair may be used
as a declared illustrative available-kit example, but is not a recovered assignment
for an actual study well. Source-era sample-sheet orientation is not a universal
rule for all current sequencers.

The RevC insert's eight-cyclei7/eight-cyclei5 programme describes the kit's standard
indexing workflow; it does not establish each study library's acquired R1 cycles or
prove its executed run recipe. PPC/NPM details, the precise full-guide revision used
by the authors, actual sample index assignments and unassembled physical end-state
modifications remain separate questions. No missing-information absence claim or
whole-protocol completion follows from this bounded positive source pass.

The source letter permits individual sequence publication with this notice:
**Oligonucleotide sequences ©2007–2012 Illumina, Inc. All rights reserved.** The whole
vendor letter is not redistributed in this repository. Retrieval URLs, external cache
names and inspected locations are recorded in `ref/MANIFEST.md`.
