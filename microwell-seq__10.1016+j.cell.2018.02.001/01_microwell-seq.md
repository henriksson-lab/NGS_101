# Microwell-seq — split-pool barcoded beads in an agarose microwell array, 3' scRNA-seq

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published, or published but not obtained. Relationships marked 🟡 *(computed)* were
> worked out with `lib/` while writing this note; they are not yet asserted in a
> self-test, because this protocol has no `tools/` module yet (status `notes`).
>
> **Source caveat.** The paper itself (Cell, Elsevier) could not be fetched, so its STAR
> Methods were **not read**. The two oligo tables (Tables S1 and S2) were read from the
> copies upstream scg_lib_structs mirrors under their Elsevier file names; their sequences
> are marked 🟢 as the paper's own supplement. Everything about reaction order and
> conditions comes from the upstream page (🟡, secondary) or is missing (🔴).

**Microwell-seq** — Han X, Wang R, Zhou Y, Fei L, Sun H, Lai S, … Orkin SH, Yuan G-C,
Chen M, Guo G. "Mapping the Mouse Cell Atlas by Microwell-Seq." *Cell* 2018;172:1091–1107.
doi:[10.1016/j.cell.2018.02.001](https://doi.org/10.1016/j.cell.2018.02.001), PMID
29474909. A published correction may exist; its DOI was not found in any fetched source and
is not given here 🔴.
Data: GEO GSE108097 (as cited by the same lab's adrenal paper, below 🟡).

Other papers listed for this method in `catalogue/scg_lib_structs.tsv`:

- nanoCAGE / CAGEscan — Plessy C et al., *Nat Methods* 2010;7:528–534,
  doi:10.1038/nmeth.1470 (PMC2906222). Cited upstream for **semi-suppressive PCR**.

Same-lab papers found while looking for open Methods (not in the catalogue):

- Lai S et al., *Cell Discov* 2018;4:34, doi:10.1038/s41421-018-0038-x — a letter that
  uses Microwell-seq; no library methods.
- Lai S et al., *Cell Regen* 2020;9:11, doi:10.1186/s13619-020-00042-8 — adrenal atlas;
  one paragraph on library construction (the transposase, below).
- Chen H et al., *Cell Discov* 2021;7:107, doi:10.1038/s41421-021-00333-7 —
  **Microwell-seq 2.0**, a different chemistry (fixed cells, well-specific RT barcodes as
  round 1). Not covered here.

Sources read (in `_data/sources/microwell-seq__10.1016+j.cell.2018.02.001/`, never
committed):

| File | What | Used for |
|---|---|---|
| `upstream_1-s2.0-S0092867418301168-mmc2.xlsx` | **Table S2**, "Oligonucleotide sequences used in library construction" (via upstream mirror) | **every library oligo** |
| `upstream_1-s2.0-S0092867418301168-mmc1.xlsx` | **Table S1**, the 3 × 96 split-pool oligos A-1…C-96 (via upstream mirror) | bead synthesis, barcode orientation |
| `upstream_Microwell_A1-A96.txt`, `_B1-B96.txt`, `_C1-C96.txt` | upstream's plain-text copies of Table S1 | checked identical to Table S1 🟡 *(computed)* |
| `upstream_Microwell-seq.html` | scg_lib_structs page | reaction order, library, read layout (🟡) |
| `same-lab_PMC7396412.xml` | adrenal paper, Europe PMC full text | transposase identity, "agarose" device |
| `same-lab_PMC6028447.xml`, `same-lab_PMC8575926.xml` | Cell Discov letter; Microwell-seq 2.0 | nothing chemical for 1.0 |
| `nmeth.1470_PMC2906222.xml` | nanoCAGE paper | what "semi-suppressive PCR" means |

Not obtained 🔴 (Cloudflare / Elsevier block from this host; fetch by hand):

| URL | Why it matters |
|---|---|
| https://www.cell.com/cell/fulltext/S0092-8674(18)30116-8 | STAR Methods: array fabrication, bead synthesis, lysis, RT, exonuclease, PCR cycles, tagmentation, sequencing read lengths |
| http://www.cell.com/article/S0092867418301168/pdf | same, as PDF |
| https://www.sciencedirect.com/science/article/pii/S0092867418301168 | supplements mmc1–mmc7 from the publisher (the two oligo tables were taken from the upstream mirror instead) |

---

## 1. What it is

A **bead-in-microwell** 3' scRNA-seq method: cells settle by gravity into an **agarose**
microwell array (agarose: 🟢 adrenal paper, discussion), one barcoded **magnetic** bead
per well is added (🟡 upstream), cells are lysed, poly(A) mRNA hybridises to the bead's
oligo-dT, and the beads are pooled for RT, template switching and whole-transcriptome
amplification. Upstream calls it "a well-based inDrop" (🟡).

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | **Cheap array**: agarose wells and ordinary magnetic beads instead of a microfluidic chip | Seq-Well (PDMS wells + Drop-seq beads), BD Rhapsody |
| 2 | Bead barcodes built by **3-round split-pool primer extension** (not ligation, not phosphoramidite split-pool) from 3 × 96 oligos | inDrop (2 rounds of extension); SPLiT-seq / sci-RNA-seq do split-pool on cells instead |
| 3 | **One 96-member 6-mer set reused in all three rounds** (96³ ≈ 8.8 × 10⁵ combinations) | — |
| 4 | Bead-end vs TSO-end discrimination by a **3'-terminal mismatch on a phosphorothioate-protected P5 primer** (§3, 🟡 inferred) | — |

The cDNA side is ordinary Drop-seq-style chemistry: SMART handle on the bead, oligo-dT
priming ([reverse transcription](../ref/concepts/reverse-transcription.md)),
[template switching](../ref/concepts/template-switching.md), single-primer PCR, then
[Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md) and a 3'-end-selecting PCR.

## 2. Oligos

### 2a. Library oligos — 🟢 verbatim from Table S2

Modification notation as the table writes it: `/rG/` = ribo-G, `/iXNA_G/` = a modified
G (the oligo is named "TSO LNA", so presumably LNA-G 🟡), `-s-` = phosphorothioate bond,
`jjjjjj` = 6-nt split-pool barcode, `nnnnnn` = 6-nt UMI.

```
indexed bead seqA  5'-Bead-Linker-TTTAGGGATAACAGGGTAATAAGCAGTGGTATCAACGCAGAGTACGTjjjjjjCGACTCACTACAGGGjjjjjjTCGGTGACACGATCGjjjjjjnnnnnnTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
indexed bead seqB  5'-Bead-Linker-TTTAGGGATAACAGGGTAATAAGCAGTGGTATCAACGCAGAGTACjjjjjjCGACTCACTACAGGGjjjjjjTCGGTGACACGATCGjjjjjjnnnnnnTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
TSO LNA            AAGCAGTGGTATCAACGCAGAGTGAAT/rG/rG/iXNA_G
TSO-PCR            AAGCAGTGGTATCAACGCAGAGT
P5                 AATGATACGGCGACCACCGAGATCTACACGCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGT-s-A-s-C
Read1 SeqA-GT      GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTACGT
Read1 SeqB -TAC    GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC
```

i7 index primers, all `P7 · <8 nt> · s7`, 47 nt (🟢 Table S2):

```
N701 CAAGCAGAAGACGGCATACGAGATTAAGGCGAGTCTCGTGGGCTCGG
N702 CAAGCAGAAGACGGCATACGAGATCGTACTAGGTCTCGTGGGCTCGG
N703 CAAGCAGAAGACGGCATACGAGATAGGCAGAAGTCTCGTGGGCTCGG
N704 CAAGCAGAAGACGGCATACGAGATTCCTGAGCGTCTCGTGGGCTCGG
N705 CAAGCAGAAGACGGCATACGAGATGGACTCCTGTCTCGTGGGCTCGG
N706 CAAGCAGAAGACGGCATACGAGATTAGGCATGGTCTCGTGGGCTCGG
N707 CAAGCAGAAGACGGCATACGAGATCTCTCTACGTCTCGTGGGCTCGG
N708 CAAGCAGAAGACGGCATACGAGATCAGAGAGGGTCTCGTGGGCTCGG
N709 CAAGCAGAAGACGGCATACGAGATGCTACGCTGTCTCGTGGGCTCGG
N710 CAAGCAGAAGACGGCATACGAGATCGAGGCTGGTCTCGTGGGCTCGG
N711 CAAGCAGAAGACGGCATACGAGATAAGAGGCAGTCTCGTGGGCTCGG
N712 CAAGCAGAAGACGGCATACGAGATGTAGAGGAGTCTCGTGGGCTCGG
N811 CAAGCAGAAGACGGCATACGAGATGCTCATGAGTCTCGTGGGCTCGG
N812 CAAGCAGAAGACGGCATACGAGATATCTCAGGGTCTCGTGGGCTCGG
N813 CAAGCAGAAGACGGCATACGAGATACTCGCTAGTCTCGTGGGCTCGG
N814 CAAGCAGAAGACGGCATACGAGATGGAGCTACGTCTCGTGGGCTCGG
N815 CAAGCAGAAGACGGCATACGAGATGCGTAGTAGTCTCGTGGGCTCGG
N816 CAAGCAGAAGACGGCATACGAGATCGGAGCCTGTCTCGTGGGCTCGG
N817 CAAGCAGAAGACGGCATACGAGATTACGCTGCGTCTCGTGGGCTCGG
N818 CAAGCAGAAGACGGCATACGAGATATGCGCAGGTCTCGTGGGCTCGG
N819 CAAGCAGAAGACGGCATACGAGATTAGCGCTCGTCTCGTGGGCTCGG
N820 CAAGCAGAAGACGGCATACGAGATACTGAGCGGTCTCGTGGGCTCGG
N821 CAAGCAGAAGACGGCATACGAGATCCTAAGACGTCTCGTGGGCTCGG
N822 CAAGCAGAAGACGGCATACGAGATCGATCAGTGTCTCGTGGGCTCGG
N823 CAAGCAGAAGACGGCATACGAGATTGCAGCTAGTCTCGTGGGCTCGG
N824 CAAGCAGAAGACGGCATACGAGATTCGACGTCGTCTCGTGGGCTCGG
```

The table also lists N801–N810; they are **byte-identical** to N701–N707, N710, N711,
N712 respectively (N708 and N709 have no N8 twin), so the table holds 36 rows but only
**26 distinct indices**. 🟡 *(computed)*

Sequencing primers given only by upstream (not in Table S2) 🟡:

```
Illumina Nextera Index 1 sequencing primer  CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
Illumina Nextera Read 2 sequencing primer   GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
```

### 2b. Bead-synthesis oligos — 🟢 Table S1 (96 of each; first member shown)

```
A-1   TTTAGGGATAACAGGGTAATAAGCAGTGGTATCAACGCAGAGTACGTTTTAGGCGACTCACTACAGGG     68 nt
B-1   CGATCGTGTCACCGACCTAAACCCTGTAGTGAGTCG                                     36 nt
C-1   AAAAAAAAAAAAAAAAAAAAAAAAAAAAAANNNNNNCCTAAACGATCGTGTCACCGA                57 nt
```

General form (🟡, from the 288 rows): `A = <bead-head> · bc1 · L1`,
`B = L2' · <bc2 as written> · L1'`, `C = A30 · N6 · <bc3 as written> · L2'`, where
`L1 = CGACTCACTACAGGG`, `L2 = TCGGTGACACGATCG` and `'` is the reverse complement.

### 2c. How they interlock — 🟡 *(computed)*

**Bead head.** `TTTAGGGATAACAGGGTAAT` (20 nt) then `rt.SMART_HANDLE` (23 nt) then `AC`
(seqB) or `ACGT` (seqA). The 20-nt head contains the 18-bp **I-SceI** homing-endonuclease
site `TAGGGATAACAGGGTAAT` (offset 2) — a rare-cutter that would release the oligo from the
bead; whether it was ever used is not said 🔴. `SMART_HANDLE + AC` is the Drop-seq bead
head; seqA adds `GT`. The two bead batches differ only by that `GT` (upstream says "two
bases" 🟡, confirmed).

**Which bead batch.** Every A-oligo in Table S1 contains `…GAGTACGT` before its barcode, so
the A-oligos build **seqA** beads, as upstream concluded. 🟡 *(computed)*

**Split-pool by primer extension.**

- B's 3' 15 nt are exactly `revcomp(L1)`, so B anneals to the 3' end of the bead-bound A
  and the polymerase copies B onto the bead strand. C's 3' 15 nt are `revcomp(L2)`, so C
  anneals to the new 3' end and copies `bc3 · UMI · T30` onto the bead.
- Tm by `chemdraw.tm` (0.5 µM, 50 mM Na⁺): L1 47 °C, L2 51 °C — 15-bp bridges.
- **Barcode orientation.** The 96 6-mers of B and C are written as the **reverse
  complement** of A's, in the same order, and B's list equals C's list exactly (e.g. A-1
  `TTTAGG`, B-1 and C-1 `CCTAAA`). After extension the bead strand therefore carries the
  **A-list sequence at all three barcode positions**, and the UMI on the bead is the
  complement of C's random `NNNNNN`. One set of 96 6-mers, **minimum pairwise Hamming
  distance 3**, all 96 distinct.
- Bead oligo, seqA: 47 + 6 + 15 + 6 + 15 + 6 + 6 + 30 = **131 nt** (seqB 129) after the
  "Linker".
- The three upstream barcode text files are identical to the Table S1 rows.

**cDNA handles.** `TSO-PCR` **is** `rt.SMART_HANDLE`. `TSO LNA` = `SMART_HANDLE · GAAT ·
rG rG G*` (30 nt) — the Smart-seq2 TSO layout with `GAAT` where Smart-seq2 has `ACAT`
(from memory of the Smart-seq2 TSO; `lib/rt.py` has no TSO constant to compare against 🟡).
After RT and template switching both cDNA ends carry the SMART handle, so one primer,
TSO-PCR, amplifies the whole cDNA (single-primer, "semi-suppressive" PCR; see §3).

**P5 and the read-1 primers.**

- `P5` = **`illumina.P5`** · `GCCTGTCCGCGG` (12 nt) · **`SMART_HANDLE`** · `AC` — 66 nt.
- `Read1 SeqA-GT` = `GCCTGTCCGCGG · SMART_HANDLE · ACGT` (39 nt);
  `Read1 SeqB-TAC` = `GCCTGTCCGCGG · SMART_HANDLE · AC` (37 nt). Each ends exactly where
  bc1 begins on its bead batch, so read 1 starts at bc1 either way. Tm 73 / 72 °C vs 59 °C
  for the bare handle — the 12-nt `GCCTGTCCGCGG` spacer, which P5 writes into the library,
  is what makes a long, high-Tm custom read-1 primer possible. Its origin is not given 🔴;
  it matches no constant in `lib/`.
- **Bead end vs TSO end.** P5's 3'-terminal 25 nt (`SMART_HANDLE + AC`) match the bead end
  of the cDNA exactly, but the TSO end reads `SMART_HANDLE + GA…`: **the last two bases of
  P5 mismatch** (A/G, C/A). Those two 3' bases are exactly the ones Table S2 protects with
  phosphorothioates, which stops a proofreading polymerase from trimming the mismatch
  away. So P5 primes only on the bead (3', barcoded) end, not on the TSO (5') end. The
  `GAAT` in the TSO (instead of an `AC`-starting tail) looks designed for this. 🟡
  (inferred from the sequences; the methods were not read)

**i7 primers.** All 36 are `illumina.P7 · 8 nt · nextera.S7` (47 nt): the standard
Nextera N7xx layout, ending at s7 without any ME bases. The 8-nt index of N701 is written
`TAAGGCGA`, i.e. **as named**; the Index 1 read (primer = `nextera.INDEX1_PRIMER`,
upstream's primer is identical 🟡 *(computed)*) reads its reverse complement,
`TCGCCTTA`. Compare ASTAR-seq's Buenrostro `Ad2.1_TAAGGCGA`, which contains `TCGCCTTA` —
the opposite orientation for the same name. Whether Table S2 lists the index names'
sequence rather than the bases actually synthesised is open 🔴.

**Read 2 primer** given upstream = `nextera.READ2_PRIMER` (= `ADAPTOR_S7`). 🟡 *(computed)*

### 2d. Upstream vs Table S2 — agreements and disagreements

| Item | Table S2 / S1 (🟢) | Upstream page (🟡) | Verdict |
|---|---|---|---|
| bead seqA/seqB fixed parts, T30 | as above | same bases, `[barcode1..3][6-bp UMI](T)30` | **agree** (T run counted: 30) |
| bead attachment | `5'-Bead-Linker-` | drawn as a bar on the 5' end; "magnetic beads" | agree; linker chemistry not given 🔴 |
| TSO 3' end | `/rG/rG/iXNA_G` | `rGrG+G` | same intent; upstream normalised XNA to LNA notation |
| P5 | ends `…GAGT-s-A-s-C` | same bases, **no phosphorothioates** | upstream drops the modification |
| i7 primers | 36 rows, 26 distinct indices, index written as named | generic "Nextera (XT) N7xx `[8-bp i7 index]`" | upstream does not list them |
| Index 1 / Read 2 primers | not listed | standard Nextera primers | upstream only |
| Tn5 adapter | not listed | "not explicitly mentioned … highly likely … s7-ME" homodimer | see §3 step 7 |
| A/B/C oligos | Table S1 | text files, identical | **agree** |

## 3. Step by step

Order 🟡 from the upstream page; quantities, temperatures and cycle numbers 🔴 (STAR
Methods not obtained).

Bead synthesis (before the experiment):

1. **Round 1**: A-1…A-96 attached to magnetic beads, one oligo per well of a 96-well plate
   (upstream: "Bind"; attachment chemistry 🔴). Bead strand:
   `bead · head · SMART_HANDLE · ACGT · bc1 · L1`.
2. **Pool, split into B-1…B-96**, primer extension: B anneals by `L1'`; the bead strand is
   extended to `… · L1 · bc2 · L2`.
3. **Pool, split into C-1…C-96**, extension: C anneals by `L2'`; bead strand becomes
   `… · L2 · bc3 · UMI · T30`.
4. **Strip the complementary strand** by heat (95 °C per upstream, which asks why not NaOH
   🟡). Beads now carry single-stranded seqA oligos.

Library:

5. **Load cells** into the agarose microwell array, then beads; **lyse**; mRNA poly(A)
   anneals to the bead T30. Beads are collected (magnet) and pooled. 🟡 (upstream)
6. **RT with template switching** on the pooled beads (MMLV-type RT + TSO LNA): first
   strand `bead · … · UMI · T30 · cDNA · CCC`, then the TSO's `GGG` pairs with the
   untemplated `CCC` and the RT copies the TSO, adding `GAAT · SMART_HANDLE'` to the
   3' end. 🟡 ([template switching](../ref/concepts/template-switching.md))
7. **Whole-transcriptome PCR with TSO-PCR alone.** Both ends carry the SMART handle, so one
   primer amplifies full-length cDNA; the bead-side 20-nt head (with the I-SceI site) is
   copied at most once (first second-strand synthesis runs to the oligo's 5' end), but the
   primer anneals at the handle, so exponentially amplified products start at the handle. Upstream labels this, after
   nanoCAGE, **semi-suppressive PCR**: when both ends carry the same sequence, short
   products fold into intramolecular panhandles and are amplified poorly, while long cDNA
   amplifies (nanoCAGE Fig. 1b, read 🟢 for the principle). Purify.
8. **Tagment** the purified cDNA. The same lab later names the enzyme: a "customized
   transposase" from the **Vazyme TruePrep DNA Library Prep Kit V2 (TD513)**, "which
   carries two identical insertion sequences" (🟢 adrenal paper, 2020 — for 2018 🟡).
   With identical adapters, and with the only Tn5-side primer being an N7xx (s7) primer,
   the insertion sequence must be **s7-ME** (`nextera.ADAPTOR_S7`), as upstream inferred.
   🟡 (inferred; adapter not published 🔴)
9. **Library PCR with P5 + one N7xx.** Three fragment classes (🟡, upstream's analysis,
   re-derived):
   - **bead end** (`SMART_HANDLE · AC · barcodes · UMI · T30 · insert · ME' · s7'`):
     P5 primes the handle end, N7xx the s7 end → **the library**;
   - **internal** (`s7 · ME · insert · ME' · s7'`): same adapter both ends, panhandle
     suppression and no P5 site → lost;
   - **TSO end** (`s7 · ME · insert · CCC · ATTC · SMART_HANDLE'`): has a handle, but P5's
     3' `AC` mismatches the `GA` that follows the handle there (§2c) → not amplified, or
     poorly.
10. Size-select / sequence (🔴 conditions).

## 4. Final library — 🟡 *(assembled from the oligos; equals upstream's drawn library)*

```
5'- P5 · GCCTGTCCGCGG · SMART_HANDLE · ACGT · bc1(6) · L1(15) · bc2(6) · L2(15) · bc3(6) · UMI(6) · T30 · <cDNA, 3' end of transcript> · ME' · s7' · i7'(8) · P7' -3'
```

- 5' adapter up to the poly(T): 122 nt; 3' adapter (`ME' · s7' · i7' · P7'`) 66 nt.
  Both strings were assembled from `lib/` constants and compared character for character
  with upstream's final library: **identical**.
- seqB beads: drop the `GT`; everything else the same.
- **No i5 index**: the P5 side carries the 12-nt spacer and handle, not an index.
  Libraries are single-indexed by i7.

## 5. Read layout / sequencing

| Read | Primer | Reads | Length needed |
|---|---|---|---|
| Read 1 | custom `Read1 SeqA-GT` (or `SeqB-TAC`) 🟢 | bc1 · L1 · bc2 · L2 · bc3 · UMI, then poly(T) | ≥ 54 nt to reach the end of the UMI 🟡 *(computed)* |
| Index 1 | `nextera.INDEX1_PRIMER` (upstream) 🟡 | i7, 8 nt (the reverse complement of the index as Table S2 writes it) | 8 |
| Read 2 | `nextera.READ2_PRIMER` (upstream) 🟡 | cDNA, from the Tn5 cut towards the poly(A) | — |

Cell barcode = bc1 + bc2 + bc3 (18 nt over 54 cycles, interrupted by two 15-nt linkers);
UMI = 6 nt. The adrenal paper confirms that barcode and UMI come from read 1 and that only
read 2 is aligned (🟢, same lab). Instrument: Illumina HiSeq is mentioned for Microwell-seq
libraries in the adrenal paper (🟡); read lengths used in the 2018 paper 🔴.

## 6. Open questions

- 🔴 The STAR Methods: bead type and attachment ("Linker"), extension enzyme and
  conditions, array geometry, lysis buffer, RT enzyme and TSO concentration, exonuclease
  step (if any), PCR cycles, tagmentation input, size selection, read lengths.
- 🔴 Whether Table S2's i7 indices are the synthesised bases or the index names (§2c);
  why N801–N810 duplicate N7xx indices; whether N811–N824 are the Illumina N714–N729 set
  (not checked — no Nextera index list in `lib/`).
- 🔴 Origin of the 12-nt `GCCTGTCCGCGG` spacer and of the I-SceI site in the bead head.
- 🟡 The P5 3'-mismatch / phosphorothioate discrimination of the TSO end is inferred from
  sequence alone.
- 🟡 Tn5 adapter = s7-ME homodimer: supported by the later same-lab statement of "two
  identical insertion sequences" and by the primers, not stated for 2018.
- 🟡 `L1 = CGACTCACTACAGGG` differs from the T7 promoter tail `CGACTCACTATAGGG` by one
  base — probably a coincidence of design; not used as a promoter here.

## 7. How this note was made (tool evaluation)

`tools/get_sources.py` saved only the upstream page; it found no open copy of the paper
(Cell/ScienceDirect return a Cloudflare challenge). The two supplementary tables were
fetched by hand from upstream's `data/Microwell-seq/` mirror (linked from the upstream
page but not followed by the tool), and three same-lab papers plus the nanoCAGE paper from
Europe PMC. `tools/doctext.py` copied the Europe PMC XML through as raw XML (tags kept);
the `.xml.txt` twins were re-made by stripping tags. `tools/scrape_primers.py` found all
Table S2 and Table S1 rows, but joined the one-sequence-per-line barcode text files into a
single 5472-nt "sequence".
