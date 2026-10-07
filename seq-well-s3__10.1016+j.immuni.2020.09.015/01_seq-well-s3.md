# Seq-Well S3 — Seq-Well with a randomly primed second strand in place of reliance on template switching

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published, or published only in a file we could not obtain. Relationships marked 🟡
> *(computed)* were worked out with `lib/` while writing this note; they are not yet
> asserted in a self-test, because this protocol has no `tools/` module yet (status
> `notes`).

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
  the upstream page only for the term "semi-suppressive PCR" (single-primer WTA).

Sources read (fetched by `tools/get_sources.py`, into
`$CHEM_DATA/sources/seq-well-s3__10.1016+j.immuni.2020.09.015/`, never committed):

| File | What | Used for |
|---|---|---|
| `j.immuni.2020.09.015_PMC7562821.html.txt` | full text (PMC), incl. STAR Methods and Key Resources Table | **every oligo** (Key Resources "Oligonucleotides"), S3 conditions, WTA PCR, Nextera PCR, sequencing |
| `j.immuni.2020.09.015_PMC7562821.xml` | JATS XML of the same | checking the oligo line breaks / stray characters (`<break/>`, `<sup>∗</sup>`) |
| `upstream_SeqWell_S3.html.txt` | scg_lib_structs page | second, complete account; checked against the paper in §2 |
| `mmc1_partial_FigS1A_legend.decoded.txt` | **hand-recovered** Figure S1A legend (see below) | step list of the S3 scheme; "random octamer" wording |
| `j.immuni.2020.09.015_PMC7562821_supplementary.zip` | Europe PMC supplement bundle — **truncated** at exactly 50,000,000 bytes | salvaged mmc2, mmc4–7 (xlsx: gene tables, no chemistry) and the first ~38 MB of mmc1.pdf |

Cross-checked against a sibling protocol's sources (read only):
`_data/sources/drop-seq-seq-well__10.1016+j.cell.2015.05.002/nmeth.4179_…_MOESM237_ESM.pdf.txt`
(Seq-Well 2017 Supplementary Table 1) for the bead oligo, which the S3 paper does not print.

Not obtained:

| File | URL | What is lost |
|---|---|---|
| **mmc1.pdf** — Document S1: Figures S1–S6 **and the Seq-Well S3 Master Protocol** (54 pages, ~75 MB) | https://pmc.ncbi.nlm.nih.gov/articles/instance/7562821/bin/mmc1.pdf (PMC download gate) | the Master Protocol: RT mix and times for S3, bead numbers, cycle numbers per sample type, read lengths. Only figure pages through ~S4 inflated from the truncated zip. 🔴 |
| mmc8.pdf — Document S2 (article + supplement) | https://pmc.ncbi.nlm.nih.gov/articles/instance/7562821/bin/mmc8.pdf | same content as mmc1 + article |
| mmc3.xlsx — Table S2, cell metadata | https://pmc.ncbi.nlm.nih.gov/articles/instance/7562821/bin/mmc3.xlsx | not chemistry |
| Shalek-lab web protocol | www.shaleklab.com (named in the paper) | the "complete, updated protocol" |

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

Not in the paper's table:

- **Bead oligo** — the paper names only the beads, "Chemgenes; MACOSKO-2011-10" in the
  methods text and "Cat#MACOSKO-2011-10B" in the Key Resources Table 🟢. The
  sequence of those beads is in Seq-Well 2017 Supplementary Table 1 (sibling sources) 🟢
  there, and on the upstream page:

```
Beads-oligo-dT-seqB  |--5'- TTTTTTTAAGCAGTGGTATCAACGCAGAGTAC[12-bp cell barcode][8-bp UMI]TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT -3'
```

  The 12-bp barcode and 8-bp UMI are stated in this paper's sequencing methods 🟢.
- **N700 index oligo** — named ("N700 Index oligo", "Nextera N700 indices") 🟢, sequence
  not given. Standard form `nextera.n7xx_primer(<i7>)` =
  `CAAGCAGAAGACGGCATACGAGAT <i7> GTCTCGTGGGCTCGG` 🟡 (as on the upstream page).
- **Read 2 and i7-index sequencing primers** — not mentioned (standard kit primers). 🔴

### Upstream (scg_lib_structs) vs. the paper — agreements and disagreements

| Oligo (upstream name) | Agreement | Note |
|---|---|---|
| TSO | **agrees**, incl. `rGrGrG` 🟢 | paper prints a stray space, `…CGCAGAG TGAAT…`; upstream has none |
| Smart PCR Primer (TSO_PCR) | **agrees** 🟢 | — |
| dN-Smart Randomer (dN-SMRT) | **agrees** 🟢 | the paper calls it "S3 Randomer" / "dn-SMART oligo"; "dN-SMRT" is upstream's name |
| New-P5-SMART PCR Hybrid Oligo (P5-TSO_Hybrid) | sequence **agrees**; upstream **omits the two phosphorothioates** `∗A∗C` | same omission as on the Drop-seq page; the PS bonds protect the 3'-terminal `AC` that makes the primer bead-end-specific (below) 🟡 (inferred) |
| Custom Read 1 Primer (Read_1_Custom_SeqB) | **agrees** 🟢 | — |
| Beads-oligo-dT-seqB | not printed in the S3 paper; **agrees** with Seq-Well 2017 Table 1 (checked in the sibling note) | 🟡 for S3 (same ChemGenes catalogue number as Seq-Well 🟢) |
| Nextera N7xx, Read 2 primer, i7 index primer | not in the paper | 🟡 upstream only; standard Nextera |
| Upstream's step "(1) Add Read 1 sequencing primer (seqA)" | label says **seqA** but the drawn primer is the **seqB** sequence (the one in the paper) | internal slip on the upstream page 🟡 |
| Upstream wording "Klenow enzyme" | paper: **Klenow Fragment (3'→5' exo⁻)**, NEB M0212L 🟢 | upstream less specific |
| Upstream step 2: RT "for all cells in one reaction" after bead removal | **agrees** with STAR Methods 🟢 | — |

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
- Upstream's read 2 and i7-index primers are exactly **`nextera.READ2_PRIMER`** and
  **`nextera.INDEX1_PRIMER`**; upstream's N7xx primer is exactly `nextera.n7xx_primer`.

Consequences:

1. After S3 every bead-bound first strand (template-switched or not) has a second strand
   beginning with the handle, so the single SMART PCR primer amplifies it (WTA). 🟡
2. The **bead end** and the **TSO/randomer end** differ right after the handle (`AC` vs
   `GA`). The P5 hybrid ends `…GCAGAGT∗A∗C`, matching only the bead end; on the other end its
   last two bases are mismatched, so in the Nextera PCR only bead-end fragments extend
   from it. 🟡 (computed; upstream says the same, step 6 "Product 2")
3. With the P5 hybrid and an N7xx as the only primers, the one exponentially amplifiable
   tagmentation product per cDNA is **bead end + s7 Tn5 end** — the 3'-most fragment. 🟡
   (agrees with upstream step 6, "Product 7")
4. The randomer's reverse complement, which ends the top strand in upstream's step 5, is
   `VNNNCCNNNTCACTCTGCGTTGATACCACTGCTT` — upstream writes this correctly. (`chemdraw.revcomp`
   passes IUPAC `B` through unchanged, so `B`→`V` was done by hand.) 🟡

## 3. Step by step

Seq-Well front end (arrays, membranes, lysis, hybridisation, bead recovery) is unchanged
from Seq-Well v1: see steps 1–7 of the
[Seq-Well note](../seq-well__10.1038+nmeth.4179/01_seq-well.md). This
paper: 10,000–15,000 cells per array, preloaded with MACOSKO-2011-10 beads, sealed with a
hydroxylated polycarbonate membrane (10 nm pores), lysis, hybridisation, bead removal. 🟢

1. **Reverse transcription** (Maxima H⁻, with the TSO). Conditions are only in the
   unobtained Master Protocol 🔴; in Seq-Well v1 it was Maxima H⁻, 4 % Ficoll, 1 mM dNTPs,
   2.5 µM TSO, 30 min RT then 52 °C 90 min 🟡 (assumed unchanged). First strand:
   `bead–TTTTTTT·handle·AC·<cbc>·<umi>·T30·<cDNA, antisense>`, ending either in `CCC` + the
   copied TSO complement (template-switched) or at a random point (not). 🟡 (Fig. S1A legend
   🟢 describes both cases; it says the RT adds a 3' overhang of three C's)
2. **Exonuclease I** on the beads, to remove excess primer (NEB M0293M). 🟢 Unextended bead
   oligos have free 3' T30 and would otherwise be primed by the randomer too. 🟡
3. **Washes**: 1 × 500 µL TE-SDS (0.5 % SDS), 2 × 500 µL TE-Tween (0.01 %). 🟢
4. **Strip the RNA**: 500 µL **0.1 M NaOH**, 5 min room temperature, end-over-end rotation —
   denatures the mRNA–cDNA hybrid; the first strand stays because it is covalently on the
   bead. Remove NaOH, wash once with TE (the text says "1 M TE"). 🟢 (mechanism 🟡)
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
7. **Tagment**: 1 ng WTA at 0.2 ng/µL + 10 µL Buffer TD + 5 µL Buffer ATM (Nextera XT),
   55 °C 5 min; 5 µL Buffer NT, 5 min RT. 🟢 (Seq-Well v1 used 750 pg.)
8. **Library PCR**: + 8 µL water, 15 µL NPM, 1 µL Custom P5 hybrid oligo, 1 µL N700 index
   oligo; 95 °C 30 s; 12 × (95 °C 10 s, 55 °C 30 s, 72 °C 30 s); 72 °C 5 min. 🟢 Primer
   concentrations not stated 🔴. SPRI 0.6× and 0.8×, Qubit, Bioanalyzer/TapeStation. 🟢

### Optimisation controls (what the paper actually compared) 🟢

Beads from one PBMC array split six ways: ± second-strand synthesis, crossed with (a) normal
RT, (b) **no TSO**, (c) heat-inactivated RT. Only "no TSO, no S3" (= Seq-Well v1 without
TSO) gave no appreciable WTA product.

## 4. Final library — 🟡 (assembled from the oligos above; agrees with upstream step 7)

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

Checked: upstream's final top strand = P5 hybrid (without PS marks) + 20 N, and ends with
exactly `nextera.ME_RC` + `nextera.S7_RC` + 8 N + `illumina.P7_RC`. 🟡 (computed)

No i5 index: the P5 end is the same custom primer for every sample; samples are
demultiplexed on the N700 (i7) index only. 🟢 ("demultiplexed using Nextera N700 indices")

Note the randomer end (`handle·GA·NNNGGNNNB…`) never appears in a sequenced library: it is
at the 5'-mRNA side of the cDNA and is lost at the 3'-end selection. 🟡 — so the random
priming position is invisible in the data, except through the shorter cDNA.

## 5. Read layout / sequencing

🟢 PBMC optimisation: NextSeq 500, 75-cycle kits (NextSeq 550 v2), 2.2 pM, **custom Read 1
primer**; skin samples: NovaSeq. Read 1 carries the 12-bp cell barcode and 8-bp UMI;
read 2 the cDNA; demultiplexing on the Nextera N700 index. Exact read lengths are not in the
main text 🔴 (Master Protocol); Seq-Well v1 used 20 / 50 / 8. 🟡

| Read | Primer | Reads | Length |
|---|---|---|---|
| Read 1 | Custom Read 1 primer (P5 hybrid minus P5) 🟢 | `<cbc 12><umi 8>` then T30 | ≥ 20 🟡 |
| Index 1 (i7) | standard, = `nextera.INDEX1_PRIMER` 🟡 | N700 index | 8 🟡 |
| Read 2 | standard, = `nextera.READ2_PRIMER` 🟡 | cDNA from the Tn5 end, sense of the mRNA | 🔴 |

Read 2 reads the transcript in the **sense** orientation, starting at a Tn5 cut inside the
3' region. 🟡 (computed from the layout)

## 6. Open questions

- 🔴 The Master Protocol (in mmc1.pdf) is lost to the truncated download: RT conditions for
  S3, whether the TSO concentration changed, read lengths, cycle numbers per sample type.
- 🔴 Why the randomer tail is `NNNGGNNNB`: the fixed internal `GG` and the non-A 3' base are
  not explained. Fig. S1A's legend calls it a "random octamer"; the tail has 7 degenerate
  positions + `GG` = 9 nt, so "octamer" matches neither count. 🟡 (counted)
- 🟡 A plausible reason for the 3'-terminal `B` (C/G/T, never A): a randomer ending in `A`
  would pair on the T30 run that every bead-bound first strand carries next to the barcode,
  priming a useless handle-to-T30 product. Not stated by the authors.
- 🔴 Whether "1 M TE" in the NaOH wash is a typo (1× TE).
- 🟡 "KAPA 5X Mastermix" vs "KAPA 2x HiFi HotStart" — presumably the 2× ReadyMix.
- 🔴 Primer concentrations in the library PCR (1 µL each, stock not stated).

## 7. How this note was made (tool evaluation)

`get_sources.py` got the PMC HTML/XML and the upstream page; all eight supplements were
behind the PMC download gate, and the Europe PMC supplementary zip it saved was **cut at
exactly 50,000,000 bytes** (no central directory, `unzip` refuses it). Stream-inflating the
local headers by hand recovered mmc2 and mmc4–7 (gene tables) and ~38 of ~75 MB of mmc1.pdf
(figures only); the Figure S1 legend there uses a re-encoded font and was decoded by a
substitution table into `mmc1_partial_FigS1A_legend.decoded.txt`. `scrape_primers.py`
found every paper oligo but split the TSO (stray space) and the P5 hybrid (stray hyphen
and stray space, `∗` phosphorothioates not parsed), and `--find` reports **0 locations for the randomer**,
which is verbatim in both files — it does not match degenerate IUPAC bases. All sequences
marked 🟢 above were confirmed with `--find` (randomer by its ACGT prefix plus grep; P5
hybrid in the paper by its two halves).
