# Plate scATAC-seq and Pi-ATAC — bulk Tn5, then FACS one nucleus per well

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. Relationships marked 🟡 *(computed)* were worked out with `lib/` while
> writing this note; they are not yet asserted in a self-test, because this protocol has
> no `tools/` module yet (status `notes`).

**Plate scATAC-seq** — Chen X, Miragaia RJ, Natarajan KN, Teichmann SA. "A rapid and
robust method for single cell chromatin accessibility profiling." *Nature
Communications* 9 (2018). doi:[10.1038/s41467-018-07771-0](https://doi.org/10.1038/s41467-018-07771-0)
(PMC6297232). Data: ArrayExpress E-MTAB-6714. Code: github.com/dbrg77/plate_scATAC-seq.

**Pi-ATAC** (protein-indexed ATAC) — Chen X, Litzenburger UM, Wei Y, Schep AN, LaGory EL,
Choudhry H, Giaccia AJ, Greenleaf WJ, Chang HY. "Joint single-cell DNA accessibility and
protein epitope profiling reveals environmental regulation of epigenomic heterogeneity."
*Nature Communications* 9 (2018). doi:[10.1038/s41467-018-07115-y](https://doi.org/10.1038/s41467-018-07115-y)
(PMC6214962).

Upstream (scg_lib_structs) draws both on one page, because the chemistry is the same:
tag in bulk, sort single nuclei, fragment and index-PCR per well.

Sources read (fetched by `tools/get_sources.py`, into
`_data/sources/plate-scatac-seq-and-pi-atac-seq__10.1038+s41467-018-07771-0/`, never committed):

| File | What | Used for |
|---|---|---|
| `s41467-018-07771-0_PMC6297232.html.txt` | plate scATAC-seq paper, full text | methods (tagmentation mix, qPCR cycle choice), data/code links |
| `s41467-018-07771-0_41467_2018_7771_MOESM1_ESM.pdf.txt` | Supplementary Information incl. **Supplementary Methods** | **step-by-step protocol and every oligo** (N7xx / S5xx table) |
| `s41467-018-07771-0_41467_2018_7771_MOESM2_ESM.pdf.txt`, `MOESM7_ESM.docx.txt` | supplementary-data descriptions | confirms Supplementary Data 1–4 are QC / annotation / motif tables (not chemistry) |
| `s41467-018-07771-0_41467_2018_7771_MOESM8_ESM.pdf.txt` | reporting summary | software versions only |
| `s41467-018-07115-y_PMC6214962.html.txt` | Pi-ATAC paper, full text | **all Pi-ATAC methods** (fixation, reverse-crosslinking buffer, PCR, sequencing) |
| `s41467-018-07115-y_41467_2018_7115_MOESM1_ESM.pdf.txt` (+ MOESM2/3/5/6/12) | Pi-ATAC supplementary figures / tables | checked: no oligo sequences |
| `s41467-018-07115-y_41467_2018_7115_MOESM4,7–11_ESM.xlsx.txt` | ChromVAR / k-mer / peak tables | not chemistry |
| `upstream_plate_and_piATAC-seq.html.txt` | scg_lib_structs method page | second account; checked below |

Not available / not read:

| What | Why |
|---|---|
| Plate scATAC Supplementary Data 1–4 (`MOESM3–6_ESM.xls`) | legacy `.xls`, no `.txt` twin made; by their descriptions (MOESM2/MOESM7) QC metrics, cell annotations, marker peaks and HOMER results — not chemistry |
| Pi-ATAC step-by-step protocol, protocols.io *private* link `https://www.protocols.io/private/F59D7D2F8FD5E57A20E039E9CF7A9785` | private link; not fetched — **the Pi-ATAC barcoding-primer sequences are in neither the paper nor its supplements** 🔴 |
| Buenrostro *et al.* 2015 (*Nature* 523:486, doi:10.1038/nature14590) primer table | Pi-ATAC ref. 3, cited for its single-cell *data* preprocessing (its bulk PCR cites Buenrostro 2013, ref. 37); not fetched |

---

## 1. What it is

Both methods are **one cell per well, no combinatorial barcode**: the cell's identity is
the pair of PCR index primers in its well (i5 × i7), exactly as in a bulk dual-indexed
Nextera library. 🟢 (upstream; implied by both papers)

The trick both rely on 🟢 (stated in the plate scATAC introduction and on the upstream page): **Tn5 tags without
fragmenting**. After tagmentation the transposase stays bound to the DNA, so nuclei
tagmented in bulk stay intact and can be FACS-sorted as single nuclei. The DNA only falls
apart in the well, when Tn5 is stripped off (SDS + heat, or proteinase K). See
[Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md).

| | Plate scATAC-seq (Teichmann lab) | Pi-ATAC (Chang lab) |
|---|---|---|
| Cells | native, unfixed; optional surface stain (e.g. CD4-PE, "TagSort") | **formaldehyde-fixed (1 %, 10 min)**, permeabilised, stained for surface **and intracellular** epitopes |
| Tagmentation | THS-seq buffer (DMF), 0.01 % digitonin, Illumina Tn5 | standard ATAC (NP-40 lysis; 0.05 % Igepal in the bulk variant), Tn5 |
| Stop | TSB (20 mM EDTA) added 1:1, i.e. 10 mM final 🟡 | 40 mM EDTA |
| Sort | DAPI+ nuclei, 96/384-well plate (index sort optional) | **index sort** (protein levels recorded per well), 96-well plates |
| Well chemistry | SDS + proteinase K lysis, 65 °C Tn5 release, Tween-20 quench | **SDS-free** reverse-crosslinking buffer (Tween-20 / Igepal + proteinase K), 65 °C overnight, 80 °C 10 min |
| Index primers | **in the well before sorting** | added with the PCR mix after reverse crosslinking |
| Index set | Nextera XT N7xx × S5xx, 24 × 16 = 384 (§2) | "barcoding primers", 96 × 90 combinations — **sequences not given** 🔴 |

What is new in each 🟡:

1. **Plate scATAC**: the lysis buffer and the index primers are pre-dispensed and frozen;
   after sorting, the plate only needs heat, Tween-20, water and PCR master mix. The
   SDS that strips Tn5 is quenched by excess Tween-20 so the PCR runs without
   purification (the paper cites the SDS/proteinase-K direct-PCR trick, ref. 14).
2. **Pi-ATAC**: a **reverse-crosslinking buffer without SDS** so fixed cells can go
   straight into PCR — the old ATAC-see buffer (1 % SDS) inhibits the polymerase 🟢 —
   and per-cell protein quantification from index FACS.

## 2. Oligos

### Plate scATAC-seq — Supplementary Methods, "Oligonucleotides sequence" 🟢

Verbatim; no modifications are given (plain DNA, 5'→3').

```
N701 CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGG
N702 CAAGCAGAAGACGGCATACGAGATCTAGTACGGTCTCGTGGGCTCGG
N703 CAAGCAGAAGACGGCATACGAGATTTCTGCCTGTCTCGTGGGCTCGG
N704 CAAGCAGAAGACGGCATACGAGATGCTCAGGAGTCTCGTGGGCTCGG
N705 CAAGCAGAAGACGGCATACGAGATAGGAGTCCGTCTCGTGGGCTCGG
N706 CAAGCAGAAGACGGCATACGAGATCATGCCTAGTCTCGTGGGCTCGG
N707 CAAGCAGAAGACGGCATACGAGATGTAGAGAGGTCTCGTGGGCTCGG
N710 CAAGCAGAAGACGGCATACGAGATCAGCCTCGGTCTCGTGGGCTCGG
N711 CAAGCAGAAGACGGCATACGAGATTGCCTCTTGTCTCGTGGGCTCGG
N712 CAAGCAGAAGACGGCATACGAGATTCCTCTACGTCTCGTGGGCTCGG
N714 CAAGCAGAAGACGGCATACGAGATTCATGAGCGTCTCGTGGGCTCGG
N715 CAAGCAGAAGACGGCATACGAGATCCTGAGATGTCTCGTGGGCTCGG
N716 CAAGCAGAAGACGGCATACGAGATTAGCGAGTGTCTCGTGGGCTCGG
N718 CAAGCAGAAGACGGCATACGAGATGTAGCTCCGTCTCGTGGGCTCGG
N719 CAAGCAGAAGACGGCATACGAGATTACTACGCGTCTCGTGGGCTCGG
N720 CAAGCAGAAGACGGCATACGAGATAGGCTCCGGTCTCGTGGGCTCGG
N721 CAAGCAGAAGACGGCATACGAGATGCAGCGTAGTCTCGTGGGCTCGG
N722 CAAGCAGAAGACGGCATACGAGATCTGCGCATGTCTCGTGGGCTCGG
N723 CAAGCAGAAGACGGCATACGAGATGAGCGCTAGTCTCGTGGGCTCGG
N724 CAAGCAGAAGACGGCATACGAGATCGCTCAGTGTCTCGTGGGCTCGG
N726 CAAGCAGAAGACGGCATACGAGATGTCTTAGGGTCTCGTGGGCTCGG
N727 CAAGCAGAAGACGGCATACGAGATACTGATCGGTCTCGTGGGCTCGG
N728 CAAGCAGAAGACGGCATACGAGATTAGCTGCAGTCTCGTGGGCTCGG
N729 CAAGCAGAAGACGGCATACGAGATGACGTCGAGTCTCGTGGGCTCGG
S502 AATGATACGGCGACCACCGAGATCTACACCTCTCTATTCGTCGGCAGCGTC
S503 AATGATACGGCGACCACCGAGATCTACACTATCCTCTTCGTCGGCAGCGTC
S505 AATGATACGGCGACCACCGAGATCTACACGTAAGGAGTCGTCGGCAGCGTC
S506 AATGATACGGCGACCACCGAGATCTACACACTGCATATCGTCGGCAGCGTC
S507 AATGATACGGCGACCACCGAGATCTACACAAGGAGTATCGTCGGCAGCGTC
S508 AATGATACGGCGACCACCGAGATCTACACCTAAGCCTTCGTCGGCAGCGTC
S510 AATGATACGGCGACCACCGAGATCTACACCGTCTAATTCGTCGGCAGCGTC
S511 AATGATACGGCGACCACCGAGATCTACACTCTCTCCGTCGTCGGCAGCGTC
S513 AATGATACGGCGACCACCGAGATCTACACTCGACTAGTCGTCGGCAGCGTC
S515 AATGATACGGCGACCACCGAGATCTACACTTCTAGCTTCGTCGGCAGCGTC
S516 AATGATACGGCGACCACCGAGATCTACACCCTAGAGTTCGTCGGCAGCGTC
S517 AATGATACGGCGACCACCGAGATCTACACGCGTAAGATCGTCGGCAGCGTC
S518 AATGATACGGCGACCACCGAGATCTACACCTATTAAGTCGTCGGCAGCGTC
S520 AATGATACGGCGACCACCGAGATCTACACAAGGCTATTCGTCGGCAGCGTC
S521 AATGATACGGCGACCACCGAGATCTACACGAGCCTTATCGTCGGCAGCGTC
S522 AATGATACGGCGACCACCGAGATCTACACTTATGCGATCGTCGGCAGCGTC
```

Per well: 2 µL of a 10 µM S5xx/N7xx mix (5 µM each) 🟢.

### How they interlock — 🟡 (computed against `lib/`)

- Every **N7xx** (24 oligos, all 47 nt) = **`illumina.P7`** + 8-nt index + **`nextera.S7`**.
  Every **S5xx** (16 oligos, all 51 nt) = **`illumina.P5`** + 8-nt index + **`nextera.S5`**.
  No primer extends into the mosaic end: they are the Nextera **XT**-style index primers,
  priming on s5/s7 alone (14 / 15 nt), unlike the Buenrostro Ad1/Ad2 primers that run
  into the ME.
- 24 × 16 = **384 index pairs**: one 384-well plate per pool. All 24 i7 and all 16 i5 are
  distinct, minimum pairwise Hamming distance 4 within each set.
- Index orientation as written in the oligo: N701 carries `TCGCCTTA` (the reverse
  complement of the Illumina i7 read `TAAGGCGA`); S502 carries `CTCTCTAT` (read as such
  on the forward-strand i5 workflow, as `ATAGAGAG` on reverse-complement i5 instruments).
- The four sequencing primers on the upstream page are exactly `nextera.READ1_PRIMER`,
  `READ2_PRIMER`, `INDEX1_PRIMER` and `INDEX2_PRIMER` (scraper match).

### Check against upstream scg_lib_structs

- **Agree** 🟡 *(computed)*: the 8-nt index that upstream lists for each of the 40
  names (N701–N729, N/S502–N/S522) is, for **all 40**, exactly the 8 nt between P7/P5
  and s7/s5 in the oligo the paper prints. Upstream's generic index-primer layouts
  (`P5 [i5] s5`, `P7 [i7] s7`) match.
- **Agree**: index primers sit in the lysis buffer before sorting; 72 °C gap fill-in
  as the first PCR step.
- **Upstream generalises**: it gives the same Nextera XT primer set for Pi-ATAC. The
  Pi-ATAC paper does **not** print its primers and quotes "96 × 90 adapter combinations"
  (8640 cells) — not a 24 × 16 set, so Pi-ATAC used a different, larger barcoding-primer
  set (presumably the Buenrostro 2015 one, its ref. 3) 🟡/🔴.
- Upstream labels its products "semi-suppressive PCR" for s5/s5 and s7/s7 fragments;
  that is upstream's (and the standard Nextera) interpretation, not stated in either
  paper 🟡.

## 3. Step by step — plate scATAC-seq

🟢 from Supplementary Methods (timestamp 15-Feb-2018) unless marked.

1. **Plate prep (day before)**: per well 2 µL 2× lysis buffer (100 mM Tris-HCl pH 8.0,
   100 mM NaCl, 40 µg/mL proteinase K, 0.4 % SDS) + 2 µL S5xx/N7xx primer mix. Seal,
   −80 °C (stable for a long time per the main text).
2. **Cells**: 5k–50k cells in BSA-coated tubes (not LoBind), pellet 500 g 4 °C 5 min,
   wash twice in cold PBS.
3. **Tagment in bulk**: 50 µL = 12.5 µL 4× THS-seq TD buffer (4× = 132 mM Tris-acetate
   pH 7.8, 264 mM K-acetate, 40 mM Mg-acetate, 64 % DMF; i.e. 33 / 66 / 10 mM, 16 % DMF
   final 🟡) + 5 µL 10× digitonin (1 µL 2 % stock + 19 µL water; 0.01 % final 🟡) + 27.5 µL water + **5 µL Illumina Tn5** (Nextera kit FC-121-1030). 37 °C,
   800 rpm, 30 min. Product: genomic fragments still held together by bound Tn5, each
   cut end carrying ME + s5 or ME + s7 on the transferred strand, 9-nt gap opposite 🟡
   ([concept](../ref/concepts/tn5-tagmentation.md)).
4. **Stop**: + 50 µL TSB (10 mM Tris-HCl pH 8.0, 20 mM EDTA), ice 10 min. Add
   100–300 µL PBS/0.5 % BSA; optional DAPI.
5. **Sort** DAPI+ single nuclei into the prepared plate; spin, seal (−80 °C for weeks
   possible).
6. **Release Tn5**: 65 °C 15 min (lid 100 °C) in SDS + proteinase K — the nucleus lyses
   and Tn5 comes off, so the DNA now physically fragments 🟡.
7. **Quench SDS**: + 4 µL 10 % Tween-20, vortex; + 2 µL water; + 10 µL 2× NEBNext
   High-Fidelity master mix → 20 µL.
8. **Index PCR**: 72 °C 10 min (gap fill-in, must precede denaturation 🟡); 98 °C 5 min;
   18 × (98 °C 10 s, 63 °C 30 s, 72 °C 20 s). The cycle number came from a qPCR
   side-experiment: 8 pre-amplification cycles, then EvaGreen qPCR, giving 8 + 10 = 18.
9. **Pool** the plate (~7.2 mL), 5 vol Buffer PB, one MinElute column on vacuum, wash
   40 mL 80 % ethanol buffer, elute 3 × 12.5 µL EB.
10. **Size select**: 0.5× SPRI (remove large) then 1.2× SPRI (keep), elute 30 µL.
11. **QC** and sequence (§5).

Discrepancy 🟡: the main-text qPCR section gives the gap fill-in as **72 °C 5 min**, the
Supplementary Methods protocol as **72 °C 10 min**.

## 3b. Step by step — Pi-ATAC (differences)

🟢 from the Pi-ATAC Methods.

1. **Fix**: 1 % formaldehyde 10 min, 0.125 M glycine 5 min, RT.
2. **Permeabilise**: ATAC lysis buffer (10 mM Tris pH 7.5, 10 mM NaCl, 3 mM MgCl₂, 0.1 %
   NP-40), spin.
3. **Stain** (where used) before tagmentation: surface and intracellular / phospho
   epitopes, 30 min per antibody, RT (e.g. CD19-PE, phospho-NF-κB p65 + secondary;
   EpCAM, CD45, HIF1α). In the K562 GATA2 and 4T1/splenocyte experiments staining was
   **after** tagmentation instead.
4. **Tagment in bulk**: standard ATAC reaction, 2× reactions per 100,000 cells (8× for
   400,000), 37 °C 30 min; stop with **40 mM EDTA**; spin, resuspend in PBS.
5. **Index sort** single cells (FACS AriaII) into 96-well plates (the Discussion also
   mentions 384-well plates) pre-filled with **20 µL
   reverse-crosslinking buffer** (composition given in the bulk 4T1 Methods paragraph, which
   the single-cell sections refer back to: 50 mM Tris-HCl pH 8.0, 0.5 % Tween-20, 0.5 % Igepal
   CA-630, proteinase K "5 ng/ml" as printed 🟢 — 🟡 probably µg/mL, unverified).
6. **Reverse crosslink**: 65 °C overnight; proteinase K off at 80 °C 10 min.
7. **Index PCR**: + 25 µL 2× NEBNext HiFi master mix + barcoding primers (two wordings:
   "2.5 µl of 25 mM barcoding primer" and "5 µl of two unique primer combinations" 🟢 —
   "25 mM" is presumably 25 µM 🟡). 72 °C 5 min; 98 °C 30 s; 20 × (98 °C 10 s, 63 °C 30 s,
   72 °C 1 min).
8. **Pool**, MinElute, elute 20 µL; **6 % PAGE**, cut > 150 bp, crush-and-soak elution
   (500 mM NaCl, 1 mM EDTA, 0.5 % SDS, 55 °C overnight), Zymo ChIP DNA Clean &
   Concentrator.

## 4. Final library — 🟡 (assembled from the oligos above; same as upstream's drawing)

Top strand, 5'→3' (standard dual-indexed Nextera):

```
5'- P5 · <i5> · s5 · ME · <genomic insert> · ME' · s7' · <i7'> · P7' -3'
```

With canonical sequences filled in (computed from `illumina` / `nextera`; i5 as written
in S5xx, i7' as written in N7xx reversed-complemented):

```
AATGATACGGCGACCACCGAGATCTACAC <i5:8> TCGTCGGCAGCGTC AGATGTGTATAAGAGACAG <insert>
CTGTCTCTTATACACATCT CCGAGCCCACGAGAC <i7':8> ATCTCGTATGCCGTCTTCTGCTTG
```

136 nt of fixed sequence (incl. the two 8-nt indices) plus the insert. The ME bases
come from the Tn5 adapter, not from the primers. Upstream's final-library drawing gives
the same segment order and sequences. Only s5…s7 fragments amplify exponentially;
s5/s5 and s7/s7 fragments are suppressed 🟡 (standard Nextera; see concept note).

For Pi-ATAC the segment list is the same if Buenrostro-style Ad1/Ad2 primers were used
(they also end in P5/P7 and anneal to s5+ME / s7+ME) 🟡; the indices themselves are 🔴.

## 5. Sequencing

- **Plate scATAC** 🟢: one 384-cell pool per lane of HiSeq 2000 or one HiSeq 2500 rapid
  run; ~1 M reads/cell, ~30,000 unique; aim ≥ 100,000 reads/cell. Read lengths 🔴 not
  stated (paired-end implied by the cutadapt trimming of Nextera sequence at the 3' end
  of short inserts 🟡).
- **Pi-ATAC** 🟢: **2 × 75** on HiSeq 4000, ~1000 cells per lane, ~0.3 M reads/cell.
- Read primers: standard Nextera — Read 1 = `nextera.READ1_PRIMER` (s5 + ME), Index 1
  (i7) = `INDEX1_PRIMER`, Index 2 (i5) = `INDEX2_PRIMER`, Read 2 = `READ2_PRIMER`
  (s7 + ME) 🟡 (upstream gives these; neither paper names read primers). Cell = (i5, i7)
  pair; no UMI, no in-read barcode.

## 6. Open questions

- 🔴 Pi-ATAC barcoding-primer sequences (96 × 90 set) — in the private protocols.io
  protocol or Buenrostro 2015; neither was read.
- 🟡 Pi-ATAC proteinase K "5 ng/ml" and primer "25 mM" are almost certainly unit typos.
- 🟡 Gap fill-in 5 min (main text) vs 10 min (Supplementary Methods) in plate scATAC.
- 🔴 Read lengths for plate scATAC.
- 🟡 Whether index FACS data were used in plate scATAC (upstream says "can be reported";
  the paper used gating — DAPI, DAPI + CD4-PE — not index values).

## 7. How this note was made (tool evaluation)

`tools/get_sources.py` fetched both papers (PMC HTML + JATS XML), all supplements twice
(Springer and Europe PMC zip copies, identical md5) and the upstream page. The four
plate-scATAC Supplementary Data files are legacy `.xls` and got no `.txt` twin (not
chemistry, so not needed). `tools/scrape_primers.py` found the upstream sequences and
the N7xx/S5xx table in the Supplementary Information, and correctly found nothing in the
Pi-ATAC files (none are printed there). A form-feed character before `S510` in the PDF
text broke a line-anchored regex in my own check script (not the tool). The methods
themselves were read, not scraped.
