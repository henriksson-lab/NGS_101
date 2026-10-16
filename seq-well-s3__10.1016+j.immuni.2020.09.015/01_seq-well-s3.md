# Seq-Well S3 — Seq-Well with a randomly primed second strand in place of reliance on template switching

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published, or published only in a file we could not obtain. Relationships marked 🟡
> *(computed)* were worked out with `lib/` and are represented by the protocol model.

**Seq-Well S3** ("Second-Strand Synthesis") — Hughes TK, Wadsworth MH II, Gierahn TM, Do T,
Weiss D, Andrade PR, Ma F, de Andrade Silva BJ, Shao S, Tsoi LC, Ordovas-Montanes J,
Gudjonsson JE, Modlin RL, Love JC, Shalek AK. "Second-Strand Synthesis-Based Massively
Parallel scRNA-Seq Reveals Cellular States and Molecular Features of Human Inflammatory
Skin Pathologies." *Immunity* 2020;53(4):878–894.e7.
doi:[10.1016/j.immuni.2020.09.015](https://doi.org/10.1016/j.immuni.2020.09.015)
(PMID 33053333, PMC7562821).

Other papers in the catalogue rows for this method:

- **Seq-Well** (the platform it modifies) — Gierahn TM et al., *Nat Methods* 2017;14:395–398,
  doi:10.1038/nmeth.4179. Written up in
  [Seq-Well](../seq-well__10.1038+nmeth.4179/01_seq-well.md).
- **Drop-seq** (same beads, same oligos) — Macosko EZ et al., *Cell* 2015, doi:10.1016/j.cell.2015.05.002.
- nanoCAGE / CAGEscan — Plessy C et al., *Nat Methods* 2010, doi:10.1038/nmeth.1470; cited by
  the original protocol for single-primer WTA.

Sources read (fetched by `tools/get_sources.py`, into
`$CHEM_DATA/sources/seq-well-s3__10.1016+j.immuni.2020.09.015/`, never committed):

| File | What | Used for |
|---|---|---|
| `j.immuni.2020.09.015_PMC7562821.html.txt` | full text (PMC), incl. STAR Methods and Key Resources Table | **every oligo** (Key Resources "Oligonucleotides"), S3 conditions, WTA PCR, Nextera PCR, sequencing |
| `j.immuni.2020.09.015_PMC7562821.xml` | JATS XML of the same | checking the oligo line breaks / stray characters (`<break/>`, `<sup>∗</sup>`) |
| `j.immuni.2020.09.015_supp_mmc1.pdf.txt` | complete Document S1, including the Seq-Well S3 Master Protocol | RT, second-strand and PCR mixes, incubation, index-primer stocks and read lengths |
| `j.immuni.2020.09.015_supp_mmc8.pdf.txt` | article plus Document S2 | cross-check of the published methods and protocol |
| `j.immuni.2020.09.015_supp_mmc3.xlsx.txt` | Table S2 cell metadata | result data, not oligo construction |
| `upstream_SeqWell_S3.html.txt` | scg_lib_structs page | orientation and historical comparison only; never evidence |
| `mmc1_partial_FigS1A_legend.decoded.txt` | **hand-recovered** Figure S1A legend (see below) | step list of the S3 scheme; "random octamer" wording |
| `j.immuni.2020.09.015_supp_mmc2.xlsx` through `mmc7.xlsx` | complete originals recovered through PMC Cloud; workbook sheets inspected | comparison statistics, cell-level metadata (including observed cell barcodes), differential expression and enrichment results; no additional synthesis recipes |

Cross-checked against a sibling protocol's sources (read only):
`/husky/henriksson/NGS101_papers/sources/drop-seq-seq-well__10.1016+j.cell.2015.05.002/nmeth.4179_…_MOESM237_ESM.pdf.txt`
(Seq-Well 2017 Supplementary Table 1) for the predecessor bead architecture.
The S3 Document S1 itself also prints SeqB in its Sequences section (extracted
lines 4718–4732); predecessor agreement supplies context rather than replacing
that direct source.

Also read the Shalek lab's [Seq-Well S³ Master Protocol](https://shaleklab.com/wp-content/uploads/2024/02/Hughes.et_.al_.Master.Protocol.pdf),
linked from its current [protocols index](https://shaleklab.com/resources/protocols/).
The 39-page author PDF is dated June 2019. Its pages 18–24 describe Exonuclease I,
NaOH stripping, second-strand synthesis and Nextera XT, and page 39 prints the
barcoded bead, modified S3 and custom Read 1 oligos. The author PDF and published
Document S1 agree on the 1xTE wash and 11 µL TD + 4 µL ATM tagmentation recipe;
the published STAR Methods prose gives different values below.

---

## 1. What it is, and what is new

Seq-Well S3 is **Seq-Well v1 with one extra step**. Everything up to and including
reverse transcription and Exonuclease I is Seq-Well (one cell + one barcoded Drop-seq bead
per sealed PDMS nanowell, lysis and mRNA capture through the membrane, beads recovered and
reverse-transcribed in bulk). Then, instead of relying solely on template switching to put
the second PCR handle on the cDNA, the beads are **stripped of RNA with 0.1 M NaOH** and a
**randomly primed second strand** is made with **Klenow exo⁻** from an oligo whose 5' part
is the SMART handle. 🟢 (Results; STAR Methods "Templated Second-Strand Synthesis")

Why: a first-strand cDNA that the RT completed but that never template-switched carries a
handle only at the bead end and cannot be amplified by the single SMART primer, so it is
lost. The random second strand puts the handle on **any** first strand, full-length or
not. 🟢 (Summary, Introduction). Consequences reported 🟢: up to 10× more transcripts and
5× more genes per cell than Seq-Well v1; S3 still yields product **without** a TSO, Seq-Well
v1 does not (Fig. S1B–E); cDNA after S3 is **shorter** than Seq-Well/Drop-seq cDNA, which
the authors list as a limitation for full-length/5'-end uses (Discussion).

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | **NaOH strip + randomly primed, handle-tagged second strand on bead-bound first-strand cDNA** (Klenow exo⁻, PEG) | the paper proposes it for Drop-seq and Slide-seq too 🟢; random-primer second-strand tagging is the general alternative to [template switching](../ref/concepts/template-switching.md) |
| 2 | Template switching kept but made **non-essential** — TSO is still in the RT; both handle sources converge on the same `handle·GA` sequence (§2) | — |

Builds on: [reverse transcription](../ref/concepts/reverse-transcription.md) from a
bead-bound oligo-dT, [template switching](../ref/concepts/template-switching.md) (SMART
handle), and [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md) (Nextera XT) with
3'-end selection by a custom P5 hybrid primer — all exactly as in
[Seq-Well v1](../seq-well__10.1038+nmeth.4179/01_seq-well.md).

## 2. Oligos

### From the paper — 🟢 Key Resources Table, "Oligonucleotides" (all "This Paper")

Written as the table gives them; `rG` = ribo-G, `∗` = phosphorothioate bond (the table uses
the U+2217 asterisk, as superscripts in the XML). The PDF/HTML wraps lines inside the
sequences and has three stray characters (a space inside the TSO, `…CGCAGAG TGAAT…`; a
hyphen `G-G` and a space `GAAGCAG TGGTA` inside the P5 hybrid; all visible in the XML as
text, not as markup); removed here.

```
Template-Switching Oligo   AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG
SMART PCR Primer           AAGCAGTGGTATCAACGCAGAGT
S3 Randomer                AAGCAGTGGTATCAACGCAGAGTGANNNGGNNNB
P5-SMART Hybrid Oligo      AATGATACGGCGACCACCGAGATCTACACGCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGT∗A∗C
Custom Read 1 Primer       GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC
```

Names used elsewhere in the same paper for the same oligos: the SMART PCR primer is the
"ISPCR oligo" in the WTA PCR paragraph, the S3 Randomer is the "dn-SMART oligo" in the
second-strand paragraph, and the P5 hybrid is the "Custom P5 hybrid Oligo" in the library
paragraph. 🟢

Not in the published article's Key Resources Table:

- **Bead oligo** — the 2019 original-author S3 Master Protocol p. 39 prints the
  SeqB bead with `JJJJJJJJJJJJNNNNNNNN` (12 cell-barcode and 8 UMI positions) and
  the flanking sequence. The published article names the beads as
  "Chemgenes; MACOSKO-2011-10" in Methods and "Cat#MACOSKO-2011-10B" in Key
  Resources. The same architecture also appears in Seq-Well 2017 Supplementary
  Table 1. 🟢 The model displays the J positions as role placeholders, because
  the source does not assign a specific bead-barcode member:

```
Beads-oligo-dT-seqB  |--5'- TTTTTTTAAGCAGTGGTATCAACGCAGAGTAC[12-bp cell barcode][8-bp UMI]TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -3'
```

  The 12-bp barcode and 8-bp UMI are also stated in this paper's sequencing methods 🟢.
- **N700 index oligo** — named ("N700 Index oligo", "Nextera N700 indices") 🟢, sequence
  of the selected member not given. The [official Illumina adapter table](https://support-docs.illumina.com/SHARE/AdapterSequences/1000000002694_22_ilmn-adapter-sequences.pdf),
  document 1000000002694 v22, Nextera DNA Indexes, prints the N7xx family; it identifies
  `TAAGGCGA` as N701, but the study does not say it used N701. The model therefore keeps
  eight variable index bases. Standard primer form `nextera.n7xx_primer(<i7>)` =
  `CAAGCAGAAGACGGCATACGAGAT <i7> GTCTCGTGGGCTCGG` 🟡 for the unverified historical kit.
- **Read 2 and i7-index sequencing primers** — model assignments from standard vendor
  Nextera context, without a historical kit-version assignment. 🟡

### Original-source boundaries

The original-author June 2019 Master Protocol p. 39 prints the complete bead, TSO,
randomer, SMART PCR, modified P5-SMART and custom Read 1 sequences. Its p. 22 names
Nextera XT, P5-SMART and N700 for the selected library. The official Illumina adapter
sequence document supplies standard Nextera family context for the right arm and
Index 1/Read 2 binding sites; the displayed bases from that context are inferred,
because the historical kit revision and chosen N700 member remain unassigned.
The upstream drawing is not evidence for these assignments.

### How the oligos interlock — 🟡 (computed with `lib/`)

One 23-nt handle, **`rt.SMART_HANDLE` = `AAGCAGTGGTATCAACGCAGAGT`**, starts or sits in
every non-Illumina oligo:

- **SMART PCR primer** *is* the handle (23 nt).
- **TSO** = handle + `GAAT` + `rGrGrG` (27 DNA + 3 RNA).
- **S3 Randomer** = handle + `GA` + `NNNGGNNNB` (34 nt). Its 3' tail is 9 nt: six `N`, a
  fixed `GG` in the middle and a 3'-terminal `B` (not A). The TSO and the randomer share
  their first **25 nt, `handle·GA`** — so whichever way the second handle got onto a cDNA,
  that end reads `handle·GA…`.
- **Bead oligo** = `TTTTTTT` + handle + `AC` + J12 + N8 + T30: the bead end reads
  `handle·AC…`.
- **P5-SMART hybrid** = **`illumina.P5`** (29) + `GCCTGTCCGCGG` (12) + handle (23) + `AC`
  (2) = 66 nt.
- **Custom Read 1 primer** = the P5 hybrid minus `illumina.P5`, exactly (37 nt); it ends on
  the bead's last constant base, so read 1 starts on cell-barcode base 1.
- Standard vendor Nextera context supplies **`nextera.READ2_PRIMER`**,
  **`nextera.INDEX1_PRIMER`** and the N7xx PCR arms; their use for the historical
  implementation is a model inference, rather than a selected index-member assignment.

Consequences:

1. After S3 every bead-bound first strand (template-switched or not) has a second strand
   beginning with the handle, so the single SMART PCR primer amplifies it (WTA). 🟡
2. The **bead end** and the **TSO/randomer end** differ right after the handle (`AC` vs
   `GA`). The P5 hybrid ends `…GCAGAGT∗A∗C`, matching only the bead end; on the other end its
   last two bases are mismatched, so in the Nextera PCR only bead-end fragments extend
   from it. 🟡 (computed from the printed bead and hybrid primer)
3. With the P5 hybrid and an N7xx as the only primers, the one exponentially amplifiable
   tagmentation product per cDNA is **bead end + s7 Tn5 end** — the 3'-most fragment. 🟡
   (conditional on the standard Nextera arm assignment)
4. The randomer's IUPAC reverse complement is
   `VNNNCCNNNTCACTCTGCGTTGATACCACTGCTT`. (`chemdraw.revcomp`
   passes IUPAC `B` through unchanged, so `B`→`V` was done by hand.) 🟡

## 3. Step by step

Seq-Well front end (arrays, membranes, lysis, hybridisation, bead recovery) is unchanged
from Seq-Well v1: see steps 1–7 of the
[Seq-Well note](../seq-well__10.1038+nmeth.4179/01_seq-well.md). This
paper: 10,000–15,000 cells per array, preloaded with MACOSKO-2011-10 beads, sealed with a
hydroxylated polycarbonate membrane (10 nm pores), lysis, hybridisation, bead removal. 🟢

1. **Reverse transcription** (Maxima H⁻, with the TSO). The Master Protocol gives a
   200 µL mix with 40 µL 5× RT buffer, 80 µL 30% PEG8000, 20 µL 10 mM dNTPs,
   5 µL 100 µM TSO, RNase inhibitor and Maxima H⁻; rotate 30 min at room temperature,
   then 90 min at 52 °C 🟢. The final TSO is 2.5 µM and PEG8000 is 12% 🟡 (computed).
   First strand:
   `bead–TTTTTTT·handle·AC·<cbc>·<umi>·T30·<cDNA, antisense>`, ending either in `CCC` + the
   copied TSO complement (template-switched) or at a random point (not). 🟡 (Fig. S1A legend
   🟢 describes both cases; it says the RT adds a 3' overhang of three C's)
2. **Exonuclease I** on the beads, to remove excess primer (NEB M0293M). 🟢 Unextended bead
   oligos have free 3' T30 and would otherwise be primed by the randomer too. 🟡 This
   occurs on the bead cDNA pool before RNA stripping, whether or not a molecule
   template-switched; the graph retains both fractions until this cleanup.
3. **Washes**: 1 × 500 µL TE-SDS (0.5 % SDS), 2 × 500 µL TE-Tween (0.01 %). 🟢
4. **Strip the RNA**: 500 µL **0.1 M NaOH**, 5 min room temperature, end-over-end rotation —
   denatures the mRNA–cDNA hybrid; the first strand stays because it is covalently on the
   bead. The published STAR Methods say to wash once with **1 M TE** (JATS line 491),
   while published Document S1 p. 19 (extracted line 3644) and the original-author
   June 2019 protocol p. 19 step 4 say **1xTE**. The graph follows the detailed
   1xTE instruction; whether the STAR Methods concentration is a transcription
   error or an intentional change remains unknown. 🟢 (mechanism 🟡)
5. **Second-strand synthesis**, 200 µL: 40 µL 5× Maxima RT buffer, 80 µL 30 % PEG 8000,
   20 µL 10 mM dNTPs, 2 µL 1 mM S3 randomer, 5 µL Klenow exo⁻, 53 µL water; **37 °C 1 h**,
   end-over-end rotation. 🟢 Final: 1× buffer, **12 % PEG**, 1 mM dNTPs, **10 µM randomer**
   🟡 (computed). The randomer's 3' tail anneals anywhere on the first strand and Klenow
   extends toward the bead, ending in the bead's handle complement; the product nearest the
   bead is always complete. Whether Klenow exo⁻ strand-displaces upstream products here is
   not discussed. 🟡 Then 2 × TE-Tween, 1 × TE, 1 × water, resuspend in 500 µL water. 🟢
6. **WTA PCR** on beads, 2,000 beads per reaction: KAPA HiFi HotStart (KK2602) 25 µL, 0.4 µL
   100 µM ISPCR (SMART PCR primer), water to 40 µL, combined with 2,000 beads 🟢. The final volume is not stated; if the
   beads come in 10 µL water (50 µL total) the primer is 0.8 µM 🟡 (assumed volume). 95 °C 3 min; 4 × (98 °C 20 s, 65 °C 45 s, 72 °C 3 min); 9–12 × (98 °C 20 s, 67 °C
   20 s, 72 °C 3 min); 72 °C 5 min. 🟢 SPRI 0.6× then 0.8×, Qubit. 🟢 (The methods call
   the mix "KAPA 5X Mastermix"; the Key Resources Table lists "KAPA 2x HiFi HotStart PCR
   mix", and Seq-Well v1 used 25 µL of 2× KAPA HiFi HotStart ReadyMix in the same 40 µL
   mastermix — so "5X" is most likely a slip for the 2× ReadyMix. 🟡)
7. **Tagment**: 1 ng WTA at 0.2 ng/µL, 55 °C 5 min; 5 µL Buffer NT, 5 min RT. 🟢
   Published STAR Methods (JATS line 509) specify 10 µL Buffer TD + 5 µL Buffer ATM;
   published Document S1 p. 22 (extracted line 3810) and the June 2019 author
   protocol p. 22 step 3 instead specify 11 µL TD + 4 µL ATM. Both total 20 µL
   with 5 µL cDNA. The graph leaves the split unassigned because the sources
   disagree; whether it changed between runs is unknown. (Seq-Well v1 used 750 pg.)
8. **Library PCR**: + 8 µL water, 15 µL NPM, 1 µL Custom P5 hybrid oligo, 1 µL N700 index
   oligo; 95 °C 30 s; 12 × (95 °C 10 s, 55 °C 30 s, 72 °C 30 s); 72 °C 5 min. 🟢 Primer
   primer stocks are both **10 µM** in the Master Protocol, giving 0.2 µM of each
   in the 50 µL final reaction 🟡 (computed). SPRI 0.6× and 0.8×, Qubit,
   Bioanalyzer/TapeStation. 🟢

### Optimisation controls (what the paper actually compared) 🟢

Beads from one PBMC array split six ways: ± second-strand synthesis, crossed with (a) normal
RT, (b) **no TSO**, (c) heat-inactivated RT. Only "no TSO, no S3" (= Seq-Well v1 without
TSO) gave no appreciable WTA product.

## 4. Final library — 🟡 (assembled from original oligos and vendor family context)

Identical in structure to Seq-Well v1; S3 changes which molecules reach it, not the layout.
Top strand 5'→3':

```
5'- P5 · GCCTGTCCGCGG · SMART handle · AC · <cbc 12> · <umi 8> · T30 · <cDNA insert> · ME' · s7' · <i7'> · P7' -3'
```

| Segment | Sequence | nt |
|---|---|---|
| `illumina.P5` | `AATGATACGGCGACCACCGAGATCTACAC` | 29 |
| spacer | `GCCTGTCCGCGG` | 12 |
| `rt.SMART_HANDLE` | `AAGCAGTGGTATCAACGCAGAGT` | 23 |
| bead constant end | `AC` | 2 |
| cell barcode | 12 | 12 |
| UMI | 8 | 8 |
| oligo-dT | T × 30 | 30 |
| cDNA insert (3' end of the mRNA, antisense on this strand) | — | variable |
| `nextera.ME_RC` | `CTGTCTCTTATACACATCT` | 19 |
| `nextera.S7_RC` | `CCGAGCCCACGAGAC` | 15 |
| i7 (reverse complement of the N7xx index) | `<i7'>` | 8 |
| `illumina.P7_RC` | `ATCTCGTATGCCGTCTTCTGCTTG` | 24 |

The P5 end is derived from the original modified hybrid primer. The Nextera right
arm is conditional on the official standard family assignment; its bases are marked
inferred in the molecular constructors and public page.

In the displayed MiSeq/NextSeq implementation there is no i5 index: the P5 end is the same custom primer for every sample; samples are
demultiplexed on the N700 (i7) index only. 🟢 ("demultiplexed using Nextera N700 indices")

Note the randomer end (`handle·GA·NNNGGNNNB…`) never appears in a sequenced library: it is
at the 5'-mRNA side of the cDNA and is lost at the 3'-end selection. 🟡 — so the random
priming position is invisible in the data, except through the shorter cDNA.

## 5. Read layout / sequencing

🟢 PBMC optimisation: NextSeq 500, 75-cycle kits (NextSeq 550 v2), 2.2 pM, **custom Read 1
primer**; skin samples: NovaSeq. Read 1 carries the 12-bp cell barcode and 8-bp UMI;
read 2 the cDNA; demultiplexing on the Nextera N700 index. The Master Protocol specifies
Read 1 **20 bp**, Read 2 **50 bp** on MiSeq/NextSeq or **50–80 bp** on NovaSeq, and an
8-bp Index 1; optional 8-bp Index 2 on NovaSeq 🟢. It cautions that Read 1 can require
21 cycles for some bead lots. 🟢

| Read | Primer | Reads | Length |
|---|---|---|---|
| Read 1 | Custom Read 1 primer (P5 hybrid minus P5) 🟢 | `<cbc 12><umi 8>` then T30 | 20; sometimes 21 🟢 |
| Index 1 (i7) | standard, = `nextera.INDEX1_PRIMER` 🟡 | N700 index | 8 🟢 |
| Read 2 | standard, = `nextera.READ2_PRIMER` 🟡 | cDNA from the Tn5 end, sense of the mRNA | 50–80 🟢 |

Read 2 reads the transcript in the **sense** orientation, starting at a Tn5 cut inside the
3' region. 🟡 (computed from the layout)

## 6. Open questions

- 🟡 Unresolved design rationale for why the randomer tail is `NNNGGNNNB`: the fixed internal `GG` and the non-A 3' base are
  not assigned by this model. Fig. S1A's legend calls it a "random octamer"; the tail has 7 degenerate
  positions + `GG` = 9 nt, so "octamer" matches neither count. 🟡 (counted)
- 🟡 A plausible reason for the 3'-terminal `B` (C/G/T, never A): a randomer ending in `A`
  would pair on the T30 run that every bead-bound first strand carries next to the barcode,
  priming a useless handle-to-T30 product. Not stated by the authors.
- 🟡 Unresolved interpretation of why the published STAR Methods say 1 M TE when both stepwise author protocols
  say 1xTE; the model follows the detailed 1xTE instruction.
- 🟡 Unresolved interpretation of why the published STAR Methods use 10 µL TD + 5 µL ATM while the 2019 and published Document S1 stepwise
  protocols use 11 µL TD + 4 µL ATM; the model has no numeric split.
- 🟡 "KAPA 5X Mastermix" vs "KAPA 2x HiFi HotStart" — presumably the 2× ReadyMix.
- 🟢 The Master Protocol supplies the library PCR primer stock concentrations: 10 µM each.

## 7. How this note was made (tool evaluation)

Initially, `get_sources.py` got the PMC HTML/XML and the upstream page. An initial Europe PMC
supplementary zip was cut at exactly 50,000,000 bytes; the complete mmc1.pdf was
subsequently recovered and its Master Protocol read. The Figure S1 legend uses a
re-encoded font and was decoded by a
substitution table into `mmc1_partial_FigS1A_legend.decoded.txt`. `scrape_primers.py`
originally split the TSO and P5 hybrid and did not match the degenerate randomer.
Those were historical tool limitations. The 2026-10-16 full-bundle scraper lists
the TSO with its RNA end, the P5 hybrid phosphorothioates and the randomer at
their original locations. The actual source table still controls transcription:
nearby labels can be joined to a hit. Complete original supplements mmc1–8 have
since been recovered through PMC Cloud. Both readers inspected every workbook's
sheet inventory and the combined article/supplement rather than relying on the
initial truncated archive or on scraper silence.

## Fresh original-source check (2026-10-16)

The retrieved originals are `Hughes.et_.al_.Master.Protocol.pdf.txt` and
`Chapter3_SeqWell.pdf.txt` in the Seq-Well S3 source bundle; the full bundle oligo
scraper was run before reading them (275 hits across ten files). The 2019 guide
lines 730–800, 930–970, 1068–1098 and 1895–1923 corroborate source order, modified
oligos, and instrument-specific read options. The 2023 chapter lines 305–330,
880–935, 1093–1170 and 1280–1310 corroborate the later 10/5 tagmentation recipe
and 20/8/50 NextSeq option. These are positive version differences, not proof that
an explanatory correction document is absent.

The 2019 guide p. 24 separately recommends optional NovaSeq Index 2 and dual
indexing, while its listed library primers are the common P5-SMART plus N700 pair.
The page therefore selects its explicit MiSeq/NextSeq 20/8/50 option; it does not
assign the optional NovaSeq second-index implementation. Bead-lot-dependent
21-cycle Read 1 is another stated variant.

Remaining linked originals need fresh source-specific checks: the 2019 guide's
NextSeq/MiSeq custom-primer guides, its generic Illumina loading-documents link,
and the 2023 chapter's cluster-density, bcl2fastq, demultiplexing and Cumulus
processing references. The chapter's scg_lib_structs illustration is excluded as
upstream artwork, not used for oligo evidence. The defining article also identifies
GEO GSE150672 and the Alexandria processed-data portal. No raw data are required
for the present molecular claims; metadata and linked-source scope remain to be
checked. This protocol remains in progress, with no new source-absence acceptance.
