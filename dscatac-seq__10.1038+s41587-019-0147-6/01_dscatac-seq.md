# dscATAC-seq / dsciATAC-seq — droplet scATAC-seq with super-loaded beads, optionally on barcoded Tn5

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. Relationships marked 🟡 *(computed)* were worked out with `lib/`; the
> schematic's construct identities are encoded in its protocol module and checked by its
> self-test. Claims that rest only on the upstream scg_lib_structs page
> (or its Bio-Rad spreadsheet) are 🟡 *(upstream)*: a secondary source, not the paper.

**dscATAC-seq / dsciATAC-seq** — Lareau CA\*, Duarte FM\*, Chew JG\*, Kartha VK, Burkett ZD,
Kohlway AS, Pokholok D, Aryee MJ, Steemers FJ, Lebofsky R, Buenrostro JD. "Droplet-based
combinatorial indexing for massive-scale single-cell chromatin accessibility."
*Nature Biotechnology* 37, 916–924 (2019). doi:[10.1038/s41587-019-0147-6](https://doi.org/10.1038/s41587-019-0147-6)
· PMID 31235917 · PMC10299900 (author manuscript). Data: GEO GSE123581.
Commercialised as the Bio-Rad **SureCell ATAC-Seq Library Prep Kit** (17004620; later
renamed ddSEQ Single-Cell ATAC) on the **ddSEQ Single-Cell Isolator**.

Other paper listed for this protocol in `catalogue/scg_lib_structs.tsv`: Plessy et al.,
"Linking promoters to functional transcripts in small samples with nanoCAGE and CAGEscan",
*Nat Methods* 2010, doi:10.1038/nmeth.1470 — cited by upstream only for the term
*semi-suppressive PCR*; not chemistry of this method, not read.

Sources read (fetched by `tools/get_sources.py`, into
`_data/sources/dscatac-seq-dsciatac-seq__10.1038+s41587-019-0147-6/`, never committed):

| File | What | Used for |
|---|---|---|
| `s41587-019-0147-6_PMC10299900.html.txt` | full text (PMC author manuscript) | **all methods**: tagmentation, droplet PCR programme, Tn5 assembly, sequencing |
| `s41587-019-0147-6_41587_2019_147_MOESM1_ESM.pdf.txt` | Supplementary Information (figures, Supplementary Note) | bead oligo cleaved off and primes at the fragment's 5' end; BAP logic |
| `s41587-019-0147-6_41587_2019_147_MOESM8_ESM.xlsx.txt` | **Supplementary Table 6** (the paper calls it Table 6) | **every published oligo**: bead structure, 96 Tn5 oligos, ME complement, N7 primers, spike-in oligo |
| `s41587-019-0147-6_41587_2019_147_MOESM7_ESM.xlsx.txt` | Supplementary Table 5 | donor / catalogue numbers only |
| `..._MOESM3..6_ESM.xlsx.txt`, `..._MOESM2_ESM.pdf.txt` | Tables 1–4 (marker genes, peaks), Reporting Summary | not chemistry |
| `upstream_dscATAC.html.txt` | scg_lib_structs page | second account: bead oligo incl. linker, kit primers, step-by-step, read layout |
| `upstream_dscATAC_dsciATAC_bead_structures_Bio-Rad_Update.xlsx.txt` | Bio-Rad bead barcode list linked from the upstream page (fetched by hand with curl) | the 192 barcodes per position and phase-block assignment |

Not obtained / not needed:

| What | URL | Why it matters |
|---|---|---|
| PMC copies of the supplements (`..._supp_NIHMS1615475-*`) | `https://pmc.ncbi.nlm.nih.gov/articles/instance/10299900/bin/NIHMS1615475-supplement-TableS4.xlsx` (and `-TableS5.tsv`, `-TableS6.tsv`; S1–S3, `-8.pdf`, `-Supplement_integrated.pdf` were saved but are reCAPTCHA pages) | redundant — the Springer `MOESM*` files are the same tables |
| SureCell ATAC-Seq Library Prep Kit User Guide (Bio-Rad 10000106678) | product page linked upstream: `http://www.bio-rad.com/en-uk/product/surecell-atac-seq-library-prep-kit?ID=PEXSR1MC1ORV` | 🔴 buffer formulations, bead oligo as manufactured, the 85 °C step, cycle table ("Table 21") — the paper defers all of these to it |

---

## 1. What it is

A **two-step, "tagment in bulk then barcode in droplets"** scATAC method 🟢: cells or
nuclei are tagmented in bulk with Tn5, kept intact, then co-encapsulated on the Bio-Rad
ddSEQ device with PCR mix and **barcoded beads**. In each droplet the bead oligos are
cleaved off the bead and act as PCR primers on the Nextera s5 end of the tagmented
fragments, putting the bead barcode at the fragment's 5' end (Supplementary Note). The
other end gets a plate-level N7xx sample index. There is **no UMI**. See
[Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md) for the adaptor chemistry.

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | **Super-loading beads**: beads loaded so densely that many droplets get ≥2; bead barcodes that share exact Tn5 insertion sites are merged afterwards (the **BAP** algorithm) 🟢 | 10x scATAC / inDrop avoid it with close-packed hydrogel beads instead |
| 2 | Bead oligo is a **primer** with a 3' s5 end, released into the droplet (by USER on a uracil linker — 🟡 *upstream*), not a capture probe | SureCell WTA (same bead platform, RNA) |
| 3 | **dsciATAC-seq**: bulk tagmentation in 24–96 wells with **barcoded Tn5** (6-nt barcode between s5 and ME), pool, then super-load *cells* too; cell = bead barcode × Tn5 barcode 🟢 | sci-ATAC-seq barcodes Tn5 the same way; SNARE/SHARE use other split-pool layers |
| 4 | Custom **Read 1 primer** so Read 1 starts directly in the bead barcode 🟢 (primer sequence 🟡 *upstream*) | Drop-seq / SureCell use a custom Read 1 too |

## 2. Oligos

### 2a. Published in the paper (Supplementary Table 6 = `MOESM8_ESM.xlsx`) — 🟢

Bead structure, as the table gives it (sheet "Bead structure"; the **"Linker sequence" cell
is empty** in the published file — checked in the raw xlsx, not a conversion loss):

```
Bead | Linker sequence | BC1 HHHHHHH | Phase Block none/A/CG/GCC/NVGC | Constant1 TATGCATGAC | BC2 HHHHHHH | Constant2 AGTCACTGAG | BC3 HHHHHHH | Nextera Read 1 TCGTCGGCAGCGTC
```

dsciATAC-seq Tn5 adaptors (sheet "dsciATAC-seq"):

```
ME complement:            /5Phos/CTGTCTCTTATACA/3ddC/
Nextera Read 2 Adapter:   GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
Nextera Read 1 Adapter (modified to contain barcode):  TCGTCGGCAGCGTC HHHHHH AGATGTGTATAAGAGACAG
  1   AAAGAA   TCGTCGGCAGCGTCAAAGAAAGATGTGTATAAGAGACAG
  2   AACAGC   TCGTCGGCAGCGTCAACAGCAGATGTGTATAAGAGACAG
  ...                                   (96 oligos, numbered 1–96, barcodes in alphabetical order)
  96  TTTGGG   TCGTCGGCAGCGTCTTTGGGAGATGTGTATAAGAGACAG
```

Sample-index primers (sheet "N70X index primer sequences"), 24 of them, N701–N707,
N710–N712, N714–N716, N718–N724, N726–N729, e.g.:

```
N701  CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGG
N729  CAAGCAGAAGACGGCATACGAGATGACGTCGAGTCTCGTGGGCTCGG
```

Spike-in for validating bead merging (sheet "Validation of multiple beads"), one
single-stranded 197-nt oligo, given in parts and as a full sequence:

```
Nextera Read 2                       GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
Constant region                      CCTAGTCGCGTAGAC
Random region                        NNNNNNNNNNNNNN
Constant region                      TACGCTAGCCCGTGGCGCGTAGCGGACTTAATCCGACCTAGCGTAGCCCGACTTAGCGTATAAGTCCGCTACGCGTCGGTCTACCCGTAGAGCGCGCTAGTCCGTGTACGGTCCGACTAG
Nextera Read 1 (reverse complement)  GACGCTGCCGACGA
```

### 2b. Only on the upstream page (Bio-Rad information to scg_lib_structs) — 🟡 *(upstream)*

```
Bead oligo   |--5'- TTTTTTTUUUTTTTTAATGATACGGCGACCACCGAGATCTACACGCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC <BC1:7> <phase:0-4> TATGCATGAC <BC2:7> AGTCACTGAG <BC3:7> TCGTCGGCAGCGTC -3'
ATAC Primer Mix   forward  AATGATACGGCGACCAC
                  reverse  CAAGCAGAAGACGGCAT
ATAC Sequencing Primer (custom Read 1)  GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC
Read 2 sequencing primer                GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
Sample Index sequencing primer          CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
Index kit (12009360)   CAAGCAGAAGACGGCATACGAGAT <8-bp sample index> GTCTCGTGGGCTCGG
```

The Bio-Rad spreadsheet linked from the upstream page gives the bead oligo **without** the
`TTTTTTTUUUTTTTT` linker (it starts `Bead-AATGATACG…`) and lists **192** 7-mers per
barcode position.

### 2c. Upstream vs paper — checked

| Item | Paper (🟢) | Upstream (🟡) | Verdict |
|---|---|---|---|
| Bead constants, phase blocks, 3' s5 | as in 2a | identical | **agree** |
| Bead linker / P5 / Read-1 handle | cell left empty | `TTTTTTTUUUTTTTT` + P5 + 37-nt custom R1 site | upstream fills a gap; not checkable against the paper 🔴 |
| Barcode alphabet "HHHHHHH" | written as H | 192 barcodes per position, **171/192 contain G** 🟡 *(computed)* | "H" is a placeholder, **not** IUPAC H (= not G). Same for the Tn5 "HHHHHH": 86/96 Tn5 barcodes contain G 🟡 *(computed)* |
| 96 Tn5 barcodes | Table 6 | Bio-Rad sheet "dsci" | **identical list, same order** 🟡 *(computed)* |
| N7 index primer | 24 explicit primers | template `P7 <8> s7` | **agree** in layout |
| Bead release | "oligonucleotides … are cleaved from the physical bead" (Suppl. Note) | USER enzyme cuts the `UUU`, 37 °C 30 min | consistent; the 37 °C 30 min step is 🟢 in the PCR programme, "USER" is not named in the paper 🟡 |
| Read lengths | R1 118, i7 8, R2 40 | same | **agree** |
| Read-1 primer | "custom sequencing primer (part of the kit)" | sequence given | upstream adds the sequence |

### How the oligos interlock — 🟡 (computed with `lib/`)

- **Tn5 adaptors.** All 96 barcoded Read-1 adaptors = **`nextera.S5` + 6-nt barcode +
  `nextera.ME`**, 39 nt each. The 96 barcodes are distinct, minimum pairwise Hamming
  distance **3** (so 1-mismatch correction is unambiguous, as the methods use).
  The Read-2 adaptor is exactly **`nextera.ADAPTOR_S7`**. The ME complement
  `CTGTCTCTTATACA` = **`nextera.ME_RC[:14]`**, and the 3' ddC sits where the 15th
  complementary base is a C — so the "15-nucleotide" ME-complementary oligo of the
  methods is 14 written bases + ddC 🟢/🟡. Mixed 1:1:2 (R1 : R2 : ME-comp), so each
  transposome half gets either the barcoded s5 arm or the plain s7 arm; only the s5 end
  carries the Tn5 barcode.
- **N7xx primers** = **`illumina.P7` + 8-nt i7 + `nextera.S7`**, 47 nt, i.e. exactly
  `nextera.n7xx_primer(i7)`. The 8 nt are written as the **reverse complement of the
  standard Illumina N7xx index** (N701 contains `TCGCCTTA`, whose reverse complement
  `TAAGGCGA` is the N701 i7 read) — true for all 24. `lib/` has no N7xx index table, so the
  comparison to Illumina's list is by name, not asserted. 🟡
- **Bead oligo** (upstream): `TTTTTTTUUUTTTTT` linker, then **`illumina.P5`** (29 nt,
  starts at position 15), then the 37-nt custom Read-1 site
  `GCCTGTCCGCGG` + **`rt.SMART_HANDLE`** + `AC` — i.e. the SMART/ISPCR handle
  `AAGCAGTGGTATCAACGCAGAGT` extended by `AC`, with a 12-nt 5' extension. Then
  BC1·phase·C1·BC2·C2·BC3 and a 3' **`nextera.S5`** (14 nt). The bead primes only through
  its 3' 14 nt, which anneal to the s5' of a tagmented fragment; a Tn5 barcode
  between s5 and ME (dsciATAC) is therefore **copied into the read, not primed over**.
- **ATAC Primer Mix**: forward = `illumina.P5[:17]`, reverse = `illumina.P7[:17]` —
  outer primers that only re-amplify complete P5…P7 products. 🟡 *(upstream sequences)*
- **Sequencing primers**: custom Read 1 = the bead's post-P5 37 nt; Read 2 primer =
  `nextera.READ2_PRIMER`; Index primer = `nextera.INDEX1_PRIMER`.
- **Barcode space** (Bio-Rad sheet): 192 7-mers per position, the **same list at BC1,
  BC2 and BC3**, minimum Hamming distance 3; the phase block is tied to the BC1 entry
  (39 none, 38 `A`, 39 `CG`, 38 `GCC`, 38 `NVGC`). 192³ ≈ 7.1 × 10⁶ bead barcodes; × 96
  Tn5 barcodes for dsciATAC. 🟡 *(computed from upstream file)*
- **Spike-in oligo**: 197 nt = `nextera.ADAPTOR_S7` + 15-nt constant + N₁₄ + 120-nt
  constant + `nextera.S5_RC`. Its 3' end is complementary to the bead oligo's 3' s5, so
  the bead primer extends on it and N7xx primes on its S7 end — it amplifies like an
  ATAC fragment, with no ME on the s5 side. Read 2 (S7+ME primer) reads the 15-nt
  constant then the 14 random bases: 29 nt, inside the 40-cycle Read 2 — matching the
  methods ("15-bp constant sequence … 14 bases downstream"). 🟡

## 3. Step by step

🟢 from the methods unless marked; formulations are deferred by the paper to the Bio-Rad
User Guide (🔴 not obtained).

1. **Cells / nuclei**: cell lines lysed Omni-ATAC style (0.1 % NP-40, 0.1 % Tween-20,
   0.01 % digitonin, 10 mM NaCl, 3 mM MgCl₂, 10 mM Tris pH 7.4, 3 min on ice), diluted
   with ATAC-Tween buffer. PBMC/BMMC: lysis during tagmentation (Whole Cell Tagmentation
   Mix, 0.1 % Tween-20, 0.01 % digitonin, PBS + 0.1 % BSA). Mouse brain: nuclei washed in ATAC-Tween buffer, then the whole-cell protocol.
2. **Bulk tagmentation**: ATAC Tagmentation Buffer + ATAC Tagmentation Enzyme (a
   concentrated Tn5 that the paper benchmarks against Illumina TDE1; 2.5 µL in 50 µL), 37 °C 30 min with agitation
   (ThermoMixer); then on ice. Products: s5/s5, s7/s7 and s5/s7 fragments with the 9-bp
   gap ([Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md)); Tn5 stays bound so
   nuclei stay intact.
   - **dsciATAC variant**: Tn5 assembled in-house — barcoded R1 adaptor : R2 adaptor :
     ME-complement 1:1:2, 100 µM total in 50 mM NaCl, 85 °C then −1 °C/min to 20 °C,
     hold 2 min; diluted 1:1 in glycerol, mixed 1:1 with 14.8 µM Tn5, 30 min RT, stored
     −20 °C. One barcode per well: 24-plex (oligos 1–3, 13–15, …, 85–87), 48-plex
     (1–6, 13–18, …, 85–90), or 96-plex (1–96; 20 µL, 4 µL Tn5, 8,000 BMMCs per well).
     Wells pooled, washed in buffer without Tn5.
3. **Droplets**: tagmented cells/nuclei + PCR mix + beads on the ddSEQ Single-Cell
   Isolator; beads super-loaded (200–5,000 beads/µL tested). Index kit 12009360 supplies
   one N7xx per sample.
4. **In-droplet programme** (C1000 Touch, deep-well module):
   37 °C 30 min → 85 °C 10 min → 72 °C 5 min → 98 °C 30 s → 8 × (98 °C 10 s, 55 °C 30 s,
   72 °C 60 s) → 72 °C 5 min. 🟢
   Interpretation 🟡 *(upstream + computed)*: 37 °C = USER cleaves the `UUU` linker and
   releases the bead oligo into the droplet; 72 °C 5 min = gap fill-in of the 9-bp
   Tn5 gaps, which also completes the s5'/s7' complements the primers anneal to;
   cycling = bead oligo (3' s5) × N7xx (3' s7) amplify s5/s7 fragments. s5/s5 fragments
   get the bead primer at both ends (P5 only) and s7/s7 the N7 primer at both ends
   (P7 only), so neither forms a clusterable library. 85 °C 10 min: purpose not stated 🔴
   (USER inactivation and/or Tn5 release are plausible 🟡).
5. **Break emulsion**, AMPure XP clean-up.
6. **Second PCR**: 98 °C 30 s → 7–9 cycles (by cell input) of 98 °C 10 s, 55 °C 30 s,
   72 °C 60 s → 72 °C 5 min. Primers not named in the paper; upstream: the kit's ATAC
   Primer Mix (P5[:17] / P7[:17]) 🟡. AMPure XP, Bioanalyzer HS.
7. **Spike-in experiment only**: the 197-nt random oligo library added at 5 nM after
   tagmentation, before encapsulation.

## 4. Final libraries — 🟡 (assembled from the oligos above)

**dscATAC-seq** (top strand, 5'→3'):

```
5'- P5 · R1-custom(37) · BC1(7) · phase(0-4) · TATGCATGAC · BC2(7) · AGTCACTGAG · BC3(7) · s5 · ME · <genomic insert> · ME' · s7' · i7'(8) · P7' -3'
```

**dsciATAC-seq**: as above with the 6-nt Tn5 barcode between s5 and ME:

```
5'- P5 · R1-custom(37) · BC1 · phase · C1 · BC2 · C2 · BC3 · s5 · <tn5bc:6> · ME · <genomic insert> · ME' · s7' · i7' · P7' -3'
```

The upstream drawing agrees with both segment lists (its top strand starts with 5
leftover T's of the linker after the first PCR, removed in practice by the second PCR with
the 17-nt P5 primer — upstream draws the final library without them). 🟡

## 5. Sequencing

🟢 NextSeq 550, High Output 150-cycle kit, loaded at 1.5 pM. **Read 1 118 cycles
(custom primer from the kit), i7 8 cycles, Read 2 40 cycles.** No i5 read.

| Read | Primer | Reads |
|---|---|---|
| Read 1 (118) | ATAC Sequencing Primer 🟡 *(upstream sequence)* | BC1 · phase · C1 · BC2 · C2 · BC3 · s5 · ME = **74–78 nt** (dsci: + Tn5 bc = **80–84 nt**), then **40–44 nt** of genome (dsci: **34–38 nt**) 🟡 *(computed)* |
| Index 1 (8) | `nextera.INDEX1_PRIMER` 🟡 | sample index (N7xx) |
| Read 2 (40) | `nextera.READ2_PRIMER` 🟡 | genome from the s7 end |

Barcodes are parsed with UMI-tools and assigned to the closest whitelist 6-/7-mer with ≤1
mismatch 🟢 — consistent with the minimum Hamming distance 3 of both lists 🟡.

## 6. Open questions

- 🔴 The bead oligo **5' linker** and **custom Read-1 site** are not in the paper (empty
  "Linker sequence" cell); they rest on Bio-Rad information relayed by upstream. The
  Bio-Rad spreadsheet itself omits the `TTTTTTTUUUTTTTT` linker.
- 🔴 Purpose of the **85 °C 10 min** step; whether USER is the release mechanism (paper
  says only "cleaved").
- 🔴 Second-PCR primers and cycle table are in the User Guide (not obtained).
- 🟡 "HHHHHHH"/"HHHHHH" in the tables is not IUPAC H — most barcodes contain G.
- 🟡 Whether the 24 N7xx primers of Table 6 are identical to the ddSEQ Index Kit
  (12009360) content: plausible, not stated.

## 7. How this note was made (tool evaluation)

`get_sources.py` got the PMC full text and all eight Springer supplements; the PMC
copies of the supplements came back as reCAPTCHA HTML pages saved under `.pdf`/`.xlsx`
names (and S4–S6 as `(manual)`), but the Springer files made them redundant. It did not
follow the upstream page's link to the Bio-Rad bead spreadsheet; that was fetched by hand
and converted with `doctext.py`. `scrape_primers.py` put Supplementary Table 6 first and
recognised P7/S7/ME in the index primers and adaptors, but did not flag `nextera.S5` in
the 96 barcoded Tn5 oligos (only `nextera.ME`). Reaction order came from reading the
methods.
