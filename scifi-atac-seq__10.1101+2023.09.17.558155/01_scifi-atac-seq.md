# scifi-ATAC-seq — pre-indexed Tn5 plus an overloaded 10x Chromium scATAC run

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. Relationships marked 🟡 *(computed)* were worked out with `lib/` while
> writing this note; they are not yet asserted in a self-test, because this protocol has
> no `tools/` module yet (status `notes`).

**scifi-ATAC-seq** — Zhang X, Marand AP, Yan H, Schmitz RJ. "scifi-ATAC-seq:
massive-scale single-cell chromatin accessibility sequencing using combinatorial fluidic
indexing." *Genome Biology* 25:90 (2024). doi:[10.1186/s13059-024-03235-5](https://doi.org/10.1186/s13059-024-03235-5)
(PMC11003106). Preprint: "Massive-scale single-cell chromatin accessibility sequencing
using combinatorial fluidic indexing", *bioRxiv* 2023,
doi:[10.1101/2023.09.17.558155](https://doi.org/10.1101/2023.09.17.558155) (v3, Feb 2024).
Data: SRA PRJNA996051; code github.com/schmitzlab/scifi-ATAC-seq.

Other papers the catalogue lists for this method: the semi-suppressive PCR reference that
upstream cites (Plessy et al., nanoCAGE/CAGEscan, *Nat Methods* 2010,
doi:10.1038/nmeth.1470). Not fetched; not needed for the chemistry.

Sources read (fetched by `tools/get_sources.py`, into
`_data/sources/scifi-atac-seq__10.1101+2023.09.17.558155/`, never committed):

| File | What | Used for |
|---|---|---|
| `s13059-024-03235-5_PMC11003106.html.txt` | Genome Biology full text | design rationale, Methods (Tn5 assembly, tagmentation, sequencing) |
| `s13059-024-03235-5_supp_..._MOESM1_ESM.pdf.txt` | Additional file 1, Supplementary Methods | **step-by-step protocol**, oligos inline |
| `s13059-024-03235-5_supp_..._MOESM2_ESM.xlsx.txt` | Additional file 2, sheet "Table S1" | **every oligo** with its 5-nt barcode (other sheets are peaks / GO tables) |
| `s13059-024-03235-5_supp_..._MOESM3_ESM.docx.txt` | peer-review history | read lengths mentioned by a reviewer; authors' clarifications |
| `2023.09.17.558155_v3.full.pdf.txt` | preprint v3 | cross-check of Methods; Fig. 1 / S1 labels |
| `2023.09.17.558155_media-1.pdf.txt` | preprint supplement | identical to MOESM1 apart from figure/table cross-reference wording 🟡 (diffed) |
| `upstream_scifi-ATAC-seq.html.txt` | scg_lib_structs page | second account of oligos, library and reads; checked below |

Not fetched / not available: the 10x Genomics user guide CG000209 (Chromium Next GEM
Single Cell ATAC v1.1, Rev E), which the paper follows "exactly" for GEM generation,
barcoding and sample-index PCR — the 10x bead oligo and the i7 sample-index primers are
therefore known here only from upstream 🟡. Get it by hand from
https://www.10xgenomics.com/support (document CG000209).

---

## 1. What it is

In short (our framing 🟡): scifi-ATAC-seq is **dsciATAC-seq-style pre-indexing on the
10x Chromium scATAC kit**; the authors say it was inspired by scifi-RNA-seq. 🟢 (main text,
paraphrased): nuclei are tagmented in a 96-well plate, each well with its own pair of
barcoded Tn5 adapters, then pooled and **overloaded** into one Chromium channel (about
100–200 k nuclei instead of the recommended maximum of 15.3 k; the Supplementary Methods
take 200–400 k into the tube, and 100 k / 200 k / 300 k inputs were compared), where the
10x bead adds the droplet barcode. A cell is the pair *(16-bp bead barcode, Tn5 well
barcode)* (Methods "Cell barcode collision detection").

| | New thing here | Builds on |
|---|---|---|
| 1 | **Two-sided** barcoded Tn5: a 5-nt barcode on the A (s5) adapter and another on the B (s7) adapter, 12 × 8 = 96 wells from 20 oligos | sci-ATAC-seq of Tu et al. 2022 (same lab); dsciATAC-seq used a one-sided barcode ([dsciATAC](../dscatac-seq__10.1038+s41587-019-0147-6/01_dscatac-seq.md)) |
| 2 | Pre-indexed nuclei loaded into the unmodified **10x scATAC v1.1** chemistry | scifi-RNA-seq ([scifi-RNA](../scifi-rna-seq__10.1101+2019.12.17.879304/01_scifi-rna-seq.md)) did the same for 10x RNA |
| 3 | A 19-nt **spacer between s5/s7 and the barcode**, so custom read and index primers are needed | — |

The chemistry of tagmentation is ordinary Tn5 ([Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md)):
the barcodes sit between the primer-entry sequence and the mosaic end, inside the
adapter, so each read starts with the Tn5 barcode then the ME.

The paper's own rationale for two-sided over one-sided barcodes 🟢 (paraphrased): (i)
fewer adapter oligos (20 for 96 combinations) and easy scaling of index complexity; (ii)
far less Tn5 (~280 µL vs >1 mL); (iii) the barcode on the s5 end helps flag and reduce
index-hopping reads (Fig. S1a,b).

## 2. Oligos

🟢 Verbatim from Additional file 2 Table S1 (sheet "Table S1"), as written there
(`[phos]` = 5' phosphate). The barcode column of the table is given after each name.

```
Tn5-ME-A1   CGTAT  5'-TCGTCGGCAGCGTCGATATGTGATAATGAGGACCGTATAGATGTGTATAAGAGACAG-3'
Tn5-ME-A2   TCGCT  5'-TCGTCGGCAGCGTCGATATGTGATAATGAGGACTCGCTAGATGTGTATAAGAGACAG-3'
Tn5-ME-A3   CTTAG  5'-TCGTCGGCAGCGTCGATATGTGATAATGAGGACCTTAGAGATGTGTATAAGAGACAG-3'
Tn5-ME-A4   CGACA  5'-TCGTCGGCAGCGTCGATATGTGATAATGAGGACCGACAAGATGTGTATAAGAGACAG-3'
Tn5-ME-A5   ACTAC  5'-TCGTCGGCAGCGTCGATATGTGATAATGAGGACACTACAGATGTGTATAAGAGACAG-3'
Tn5-ME-A6   TATAC  5'-TCGTCGGCAGCGTCGATATGTGATAATGAGGACTATACAGATGTGTATAAGAGACAG-3'
Tn5-ME-A7   GGAAT  5'-TCGTCGGCAGCGTCGATATGTGATAATGAGGACGGAATAGATGTGTATAAGAGACAG-3'
Tn5-ME-A8   AAGTC  5'-TCGTCGGCAGCGTCGATATGTGATAATGAGGACAAGTCAGATGTGTATAAGAGACAG-3'
Tn5-ME-A9   AGTAA  5'-TCGTCGGCAGCGTCGATATGTGATAATGAGGACAGTAAAGATGTGTATAAGAGACAG-3'
Tn5-ME-A10  AATTG  5'-TCGTCGGCAGCGTCGATATGTGATAATGAGGACAATTGAGATGTGTATAAGAGACAG-3'
Tn5-ME-A11  AGATG  5'-TCGTCGGCAGCGTCGATATGTGATAATGAGGACAGATGAGATGTGTATAAGAGACAG-3'
Tn5-ME-A12  CGCGA  5'-TCGTCGGCAGCGTCGATATGTGATAATGAGGACCGCGAAGATGTGTATAAGAGACAG-3'

Tn5-ME-B1   TCGGA  5'-GTCTCGTGGGCTCGGTGAATGTGTAGAAGACAGATCGGAAGATGTGTATAAGAGACAG-3'
Tn5-ME-B2   GTTTC  5'-GTCTCGTGGGCTCGGTGAATGTGTAGAAGACAGAGTTTCAGATGTGTATAAGAGACAG-3'
Tn5-ME-B3   AGCTT  5'-GTCTCGTGGGCTCGGTGAATGTGTAGAAGACAGAAGCTTAGATGTGTATAAGAGACAG-3'
Tn5-ME-B4   TCATT  5'-GTCTCGTGGGCTCGGTGAATGTGTAGAAGACAGATCATTAGATGTGTATAAGAGACAG-3'
Tn5-ME-B5   GCTCC  5'-GTCTCGTGGGCTCGGTGAATGTGTAGAAGACAGAGCTCCAGATGTGTATAAGAGACAG-3'
Tn5-ME-B6   ACTAA  5'-GTCTCGTGGGCTCGGTGAATGTGTAGAAGACAGAACTAAAGATGTGTATAAGAGACAG-3'
Tn5-ME-B7   CGAGG  5'-GTCTCGTGGGCTCGGTGAATGTGTAGAAGACAGACGAGGAGATGTGTATAAGAGACAG-3'
Tn5-ME-B8   ACCGG  5'-GTCTCGTGGGCTCGGTGAATGTGTAGAAGACAGAACCGGAGATGTGTATAAGAGACAG-3'

Tn5MErev           5'-[phos]CTGTCTCTTATACACATCT-3'

scifi-qPCR-F       5'-TCGTCGGCAGCGTCGATATGTGATAATGAGGAC-3'
scifi-qPCR-R       5'-GTCTCGTGGGCTCGGTGAATGTGTAGAAGACAGA-3'

1_Read1            5'-TCGTCGGCAGCGTCGATATGTGATAATGAGGAC-3'     HPLC purified
2_Index1(i7)       5'-TCTGTCTTCTACACATTCACCGAGCCCACGAGAC-3'    HPLC purified
3_Index2(i5)       5'-GTCCTCATTATCACATATCGACGCTGCCGACGA-3'     HPLC purified
4_Read2            5'-GTCTCGTGGGCTCGGTGAATGTGTAGAAGACAGA-3'    HPLC purified
```

The bottom (ME-complement) oligo is named three ways: "Tn5MErev" (Table S1),
"Tn5-ME-Rev" (Supplementary Methods), "Tn5-ME-rev" (main Methods). Same oligo. 🟡

Not in the paper, from upstream only (10x kit components) 🟡:

```
10x Genomics bead oligo (scATAC)        5'- AATGATACGGCGACCACCGAGATCTACAC [16-bp GEM barcode] TCGTCGGCAGCGTC -3'
i7 Sample Index Plate N, Set A          5'- CAAGCAGAAGACGGCATACGAGAT [8-bp sample index] GTCTCGTGGGCTCGG -3'
```

### How the oligos interlock — 🟡 (computed with `lib/nextera`, `lib/illumina`)

- Every A adapter = **`nextera.S5`** (14 nt) + 19-nt spacer **`GATATGTGATAATGAGGAC`** +
  5-nt barcode + **`nextera.ME`** — 57 nt. Every B adapter = **`nextera.S7`** (15 nt) +
  19-nt spacer **`TGAATGTGTAGAAGACAGA`** + 5-nt barcode + `nextera.ME` — 58 nt. The
  barcode embedded in each sequence equals the table's barcode column for all 20 oligos.
- **Tn5MErev = `nextera.ME_RC`** — the standard 5'-phosphorylated 19-nt ME complement,
  shared by both adapter types.
- **scifi-qPCR-F = 1_Read1** = S5 + spacer A (33 nt); **scifi-qPCR-R = 4_Read2** = S7 +
  spacer B (34 nt). The same oligo serves as qPCR primer and as custom sequencing primer.
  Both end immediately 5' of the barcode, so Read 1 and Read 2 each **start with the
  5-nt Tn5 barcode, then the 19-nt ME**, then genomic DNA. 
- **2_Index1(i7) = revcomp(4_Read2)** and **3_Index2(i5) = revcomp(1_Read1)**, exactly.
  The custom Index 1 primer therefore ends in `nextera.S7_RC` and reads the i7 sample
  index straight after the spacer; it replaces the standard Nextera index-1 primer
  (`ME_RC + S7_RC`), which would not anneal because the barcode and spacer sit between
  ME and s7.
- Barcode spacing: minimum Hamming distance 2 within the 12 A barcodes and within the
  8 B barcodes (so single errors are detectable but not correctable). A and B are read in
  different reads, so A-vs-B distance (min 1) does not matter.
- The 10x bead oligo ends in **`nextera.S5`**, which is exactly the first 14 nt of every A
  adapter: the bead primes from A ends only. The 10x i7 primer ends in **`nextera.S7`**,
  the first 15 nt of every B adapter: it primes from B ends only. 🟡 (upstream sequences)

## 3. Step by step

Plant (maize seedling) protocol; concentrations and times 🟢 from Supplementary Methods
unless noted.

1. **Anneal adapters**: each Tn5-ME-Ax / Tn5-ME-Bx with Tn5MErev, 1:1 at 100 µM;
   98 °C 2 min, then −1 °C per 10-s cycle to 25 °C. (Main Methods instead: 95 °C, −1 °C
   per minute to 20 °C — a disagreement, see §6.)
2. **Assemble Tn5** (in-house Tn5, Tu et al. 2020, Addgene 127916): A adapters 2 µL +
   15 µL Tn5; B adapters 3 µL + 25 µL diluted Tn5; 25 °C 1 h.
3. **Plate the combinations**: 1.5 µL A-Tn5 and 1.5 µL B-Tn5 per well of a 96-well
   plate → 3 µL Tn5 per well, a unique (A, B) pair per well.
4. **Nuclei**: chop seedlings in NIB-cutting buffer (MES-KOH pH 5.4, NaCl, sucrose,
   spermine, spermidine, DTT, BSA, 0.5 % Triton X-100), filter 40 / 20 µm, 35 % Percoll
   cushion, wash and resuspend in TAPS buffer (TAPS pH 8.0, MgCl₂) with 0.1 % Tween 20 +
   0.01 % digitonin; dilute to ~0.5–1 k nuclei/µL.
5. **Tagment**: 10 µL nuclei per well, 37 °C 60 min. Each fragment end gets an A or a B
   adapter at random, leaving the usual 9-nt gap
   ([Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md)): A–A, B–B and A–B
   fragments. 🟡
6. **Stop**: 12 µL stop buffer per well (final 10 mM Tris-HCl pH 7.8, 20 mM EDTA, 2 %
   BSA) — EDTA chelates Mg²⁺. 🟢
7. **Pool** all wells, pellet 500 rcf, resuspend in 10x Diluted Nuclei Buffer, filter
   40 µm, count, take 200–400 k nuclei, pellet, leave ~5 µL and add 7 µL 10x ATAC
   Buffer B.
8. **10x Chromium scATAC v1.1** (Chip H), "exactly" per CG000209 for GEM generation,
   barcoding and post-GEM cleanup. Inside each GEM 🟡 (from the kit as drawn by upstream,
   not stated by the paper): 72 °C 5 min fills the 9-nt gaps and displaces Tn5MErev, then
   linear PCR with the bead oligo, whose 3' `S5` anneals to the complement of the A
   adapter's s5 — adding P5 + 16-bp GEM barcode to every A end.
9. **qPCR to set cycles**: 1 µL of the cleaned GEM product in 10 µL Luna Universal qPCR
   mix with scifi-qPCR-F/R (0.5 µL each of 10 µM); 95 °C 30 s; 30 × (95 °C 15 s, 63 °C
   30 s, 72 °C 1 min).
10. **Sample-index PCR** on the remaining 40 µL, per the 10x guide, cycle number chosen
    around the Cq; aim ~100 ng, "without more than 9 additional cycles" 🟢. Primers
    (from the kit) 🟡: a P5 primer plus one i7 Sample Index primer (P7 + i7 + S7), which
    primes from B ends. Only fragments with an A end (bead-tagged, P5) and a B end (S7)
    grow exponentially; A–A and B–B fragments do not. Upstream attributes this to
    semi-suppressive PCR. 🟡

## 4. Final library — 🟡 (assembled from Table S1 plus the upstream 10x sequences)

Top strand, 5'→3' (192 nt of adapter + insert):

```
P5 (29) · GEM barcode (16) · s5 (14) · spacer A (19) · Tn5 barcode A (5) · ME (19) ·
  <genomic insert> ·
ME' (19) · Tn5 barcode B' (5) · spacer B' (19) · s7' (15) · i7' (8) · P7' (24)
```

- spacer A = `GATATGTGATAATGAGGAC`; spacer B' = revcomp(`TGAATGTGTAGAAGACAGA`).
- Barcode B appears on the top strand as the reverse complement of the Table S1 barcode;
  Read 2 reads it in the table's orientation. 🟡 (computed)

Agreement with upstream: the final structure upstream draws matches this segment for
segment, and every oligo upstream lists (Tn5-ME-A/B, bottom ME, four sequencing primers)
matches Table S1 exactly. 🟡 (checked) Upstream calls the ME complement "Tn5ME-bottom".

## 5. Sequencing

🟢 NovaSeq 6000 S4, dual-index mode, with the four custom primers (1_Read1, 2_Index1,
3_Index2, 4_Read2). Because every read starts with a 5-nt barcode plus the constant 19-nt
ME, scifi-ATAC-seq should be **< 50 % of a lane**, or spiked with PhiX. Depth: on average
7,617 read pairs per cell (62.3 % unique).

Read layout 🟡 (from the structure; upstream agrees):

| Read | Primer | Cycles | Content |
|---|---|---|---|
| Read 1 | 1_Read1 | not stated 🔴 (upstream: ≥74) | Tn5 barcode A (5) · ME (19) · gDNA |
| Index 1 (i7) | 2_Index1(i7) | 8 | 10x sample index |
| Index 2 (i5) | 3_Index2(i5) | 16 | **GEM bead barcode** |
| Read 2 | 4_Read2 | not stated 🔴 (upstream: ≥74) | Tn5 barcode B (5) · ME (19) · gDNA |

A reviewer describes reads of "R1:49bp and R2 50bp" — presumably the genomic part after
trimming; the authors did not give cycle numbers. 🟡 Processing 🟢: the 16-bp i5 bead
barcode is put in the read name with UMI-tools (`--bc-pattern=NNNNNNNNNNNNNNNN`), the
inline Tn5 barcodes are demultiplexed with cutadapt, and the cell barcode is bead + Tn5
barcode.

## 6. Open questions and disagreements between sources

- 🔴 **Rows vs columns.** Both Methods texts put the 12 A adapters "by rows" and the 8 B
  adapters "by columns", but a 96-well plate has 8 rows and 12 columns, and the Results
  say "8 rows × 12 columns … 15 µL in 8 rows and 10 µL in 12 columns". The layout that
  fits the plate is A in columns, B in rows; the text is inconsistent.
- 🔴 **Number of B adapters.** Supplementary Methods step II.1 says "12 Tn5-ME-Ax and
  18 Tn5-ME-Bx"; Table S1 and the rest of the text have 8 B adapters (B1–B8). Most likely
  a typo for 8. 🟡
- 🔴 **Tn5 assembly.** The main Methods say the annealed A and B adapters are mixed 1:1
  and 0.143 µL of the mix is added to 10 µL Tn5 — which would put A and B on the same
  transposome stock and is incompatible with the 12 × 8 combinatorial plating. The
  Supplementary Methods assemble A and B separately (2 µL + 15 µL Tn5; 3 µL + 25 µL
  diluted Tn5). Followed the supplement here.
- 🟡 **Annealing programme** differs (98 °C / −1 °C per 10 s to 25 °C vs 95 °C / −1 °C per
  min to 20 °C). Either works for a 19-bp duplex.
- 🔴 Read 1 / Read 2 cycle numbers; the exact SI-PCR primer names and cycles (deferred to
  the 10x guide, not fetched).
- 🔴 Origin of the 19-nt spacers (presumably from the lab's earlier sci-ATAC-seq design,
  Tu et al. 2022, *Plant Commun.*; not checked).
- 🟡 "Index hopping" in this paper refers to cross-contamination between nuclei sharing a
  droplet (Fig. S1a; figure not read in detail). The authors credit the s5-side barcode
  with helping to distinguish such reads, and quantify the residual contamination from
  genotype SNPs (e.g. 1.93 % in the eight-sample multiplex) 🟢; the molecular mechanism is
  not spelled out.

## 7. How this note was made (tool evaluation)

`tools/get_sources.py` had already fetched the preprint, its supplement, the PMC full
text and the Europe PMC supplementary zip. The oligo table is the first sheet of a 46 MB
xlsx whose other sheets are peak / ACR tables; its 64 MB text twin made
`tools/scrape_primers.py` (both the scan and `--find`) run for minutes without output on the
source directory. Work-around: the scan and `--find` were run on a scratch copy holding the
supplement texts, the upstream page and only the Table S1 sheet; every 🟢 sequence was
located there (and checked as an exact substring of Table S1). Oligos were taken from the "Table S1" sheet and their
relationships computed in a scratch script against `lib/nextera` and `lib/illumina`.
