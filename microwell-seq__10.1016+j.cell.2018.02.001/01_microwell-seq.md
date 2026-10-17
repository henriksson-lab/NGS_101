# Microwell-seq — split-pool barcoded beads in an agarose microwell array, 3' scRNA-seq

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published, or published but not obtained. Relationships marked 🟡 *(computed)* were
> worked out with `lib/`; the model is `tools/microwell_seq.py`.
>
> **The scg_lib_structs drawn page is not a source here.** It is inspiration only: no
> sequence, length or reaction in this note or in the model is taken from it, and
> agreement with it is not verification.
>
> **What is, and is not, published.** Every concrete base the model draws comes from the
> authors' **own Supplementary Tables S1 and S2** 🟢 — Elsevier's supplementary files
> `1-s2.0-S0092867418301168-mmc1.xlsx` and `mmc2.xlsx` for this article (PII
> S0092867418301168), mirrored in the scg_lib_structs repository under
> `data/Microwell-seq/`. A filename prefix is not provenance: those two files are the
> publisher's own supplement wherever they are read from, and `mmc2` is titled
> "Supplementary table 2  Oligonucleotide sequences used in library construction". The
> rest of the model is canonical `lib/illumina.py` / `lib/nextera.py` constants and
> placeholders. The paper's **STAR Methods were read on 2026-10-17** (hand-downloaded
> `j.cell.2018.02.001.pdf`, see §10); they supply reaction order and conditions. The
> Tn5 insertion sequence is still not printed, so the s7-ME end is an explicit inference
> drawn dotted on the page. Sections 2–6 below predate that reading; §10 records what
> the STAR Methods confirm or change, and the model follows §10.

**Microwell-seq** — Han X, Wang R, Zhou Y, Fei L, Sun H, Lai S, … Orkin SH, Yuan G-C,
Chen M, Guo G. "Mapping the Mouse Cell Atlas by Microwell-Seq." *Cell* 2018;172:1091–1107.
doi:[10.1016/j.cell.2018.02.001](https://doi.org/10.1016/j.cell.2018.02.001), PMID
29474909. Correction: *Cell* 2018;173:1307, doi:10.1016/j.cell.2018.05.012 🟢 (appended to
the PDF; fixes an author name, the bead supplier and the lysis-buffer Tris concentration,
0.1 M not 1 M).
Data: GEO GSE108097 (as cited by the same lab's adrenal paper, below 🟡).

Other papers listed for this method in `catalogue/scg_lib_structs.tsv`:

- nanoCAGE / CAGEscan — Plessy C et al., *Nat Methods* 2010;7:528–534,
  doi:10.1038/nmeth.1470 (PMC2906222). Read 🟢 for what **semi-suppressive PCR** means (its Fig. 1b); it is the origin of the single-primer design, not a source for any Microwell-seq base.

Same-lab papers found while looking for open Methods (not in the catalogue):

- Lai S et al., *Cell Discov* 2018;4:34, doi:10.1038/s41421-018-0038-x — a letter that
  uses Microwell-seq; no library methods.
- Lai S et al., *Cell Regen* 2020;9:11, doi:10.1186/s13619-020-00042-8 — adrenal atlas;
  one paragraph on library construction (the transposase, below).
- Chen H et al., *Cell Discov* 2021;7:107, doi:10.1038/s41421-021-00333-7 —
  **Microwell-seq 2.0**, a different chemistry (fixed cells, well-specific RT barcodes as
  round 1). Not covered here.

Sources read (in `/husky/henriksson/NGS101_papers/sources/microwell-seq__10.1016+j.cell.2018.02.001/`, never
committed):

| File | What | Used for |
|---|---|---|
| `1-s2.0-S0092867418301168-mmc2.xlsx` | the article's **Supplementary Table S2**, "Oligonucleotide sequences used in library construction" (Elsevier ESM; on disk from the scg_lib_structs mirror of `data/Microwell-seq/`) | **every library oligo** 🟢 |
| `1-s2.0-S0092867418301168-mmc1.xlsx` | the article's **Supplementary Table S1**, the 3 × 96 split-pool oligos A-1…C-96 (same provenance) | bead synthesis, barcode orientation 🟢 |
| `Microwell_A1-A96.txt`, `_B1-B96.txt`, `_C1-C96.txt` | plain-text splits of Table S1 kept in the same mirror | row-for-row identical to Table S1 🟡 *(computed)*; Table S1 itself is the authority |
| `upstream_Microwell-seq.html` | scg_lib_structs method page | **inspiration only, not a source** — nothing is cited from or verified against it |
| `same-lab_PMC7396412.xml` | adrenal paper, Europe PMC full text | transposase identity, "agarose" device |
| `same-lab_PMC6028447.xml`, `same-lab_PMC8575926.xml` | Cell Discov letter; Microwell-seq 2.0 | nothing chemical for 1.0 |
| `nmeth.1470_PMC2906222.xml` | nanoCAGE paper | what "semi-suppressive PCR" means |
| `lab_MCA_microwell-seq_protocol_page.html` | Guo-lab "Microwell-seq protocol" section of https://bis.zju.edu.cn/MCA/ (undated web page, retrieved 2026-10-17; §9) | lysis, RT, Exo I, WTA and library PCR conditions, Nextera XT, NextSeq 500, Table 1 primers 🟢 (lab protocol; version not stated) |
| `j.cell.2018.02.001.pdf` (+ `.txt`) | the defining *Cell* article PDF with STAR Methods and the appended correction (hand-downloaded from cell.com on the institutional network, 2026-10-17; printed DOI/title/authors checked) | Key Resources Table oligos (lines 2510–2565), Method Details (lines 2670–2740), data processing (2815–2820), correction (6335–6355) 🟢 |
| `lab_Protocol_for_making_microwells_and_barcoded_beads.docx` | Guo-lab docx linked from the same page (undated; §9) | microwell geometry, EDC coupling, Phanta extension rounds, 95 °C strand stripping 🟢 (lab protocol) |

Not obtained 🔴 (Cloudflare / Elsevier block from this host; fetch by hand):

| URL | Why it matters |
|---|---|
| https://www.sciencedirect.com/science/article/pii/S0092867418301168 | supplements mmc3–mmc7 direct from the publisher (mmc1 and mmc2, the two oligo tables, are on disk via the scg_lib_structs mirror of `data/Microwell-seq/`; `catalogue/extra_sources.tsv` now declares that mirror so `get_sources.py` fetches them as `scg_*` rather than by hand) |

---

## 1. What it is

A **bead-in-microwell** 3' scRNA-seq method: cells settle by gravity into an **agarose**
microwell array (agarose: 🟢 adrenal paper, discussion), one barcoded **magnetic** bead
per well is added (🔴 bead type and loading not read; the beads are magnetic per the
adrenal paper's retrieval step 🟢), cells are lysed, poly(A) mRNA hybridises to the bead's
oligo-dT, and the beads are pooled for RT, template switching and whole-transcriptome
amplification 🟡 *(derived from the oligo architecture: one bead-borne oligo-dT with a
SMART handle, one TSO, one TSO-PCR primer)*.

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

Table S2 lists **no** Index 1 or Read 2 sequencing primer 🔴. Because the only Tn5-side
oligo it prints is the stock N7xx layout and no custom read primer appears, the model uses
the canonical Illumina/Nextera primers `nextera.INDEX1_PRIMER` and `nextera.READ2_PRIMER`
from `lib/nextera.py` 🟡 *(derived; the bases come from Illumina's adapter document, not
from any Microwell-seq statement)*.

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
(seqB) or `ACGT` (seqA) 🟢. The 20-nt head contains the 18-bp **I-SceI** homing-endonuclease
site `TAGGGATAACAGGGTAAT` (offset 2) — a rare-cutter that would release the oligo from the
bead; whether it was ever used is not said 🔴. `SMART_HANDLE + AC` is also the Drop-seq/Seq-Well bead
head; seqA adds `GT`. The two bead designs Table S2 prints differ only by that `GT` 🟡
*(computed from the two rows)*.

**Which bead batch.** Every A-oligo in Table S1 contains `…GAGTACGT` before its barcode, so
the A-oligos build **seqA** beads. 🟡 *(computed)*

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
- The mirror's three plain-text barcode files are row-for-row identical to Table S1.

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
`lib/nextera.py`) reads its reverse complement,
`TCGCCTTA`. Compare ASTAR-seq's Buenrostro `Ad2.1_TAAGGCGA`, which contains `TCGCCTTA` —
the opposite orientation for the same name. Whether Table S2 lists the index names'
sequence rather than the bases actually synthesised is open 🔴.

**Read 2 primer**: `nextera.READ2_PRIMER` (= `ADAPTOR_S7`), derived as in §2c 🟡; Table S2
prints no read primer for this end 🔴.

### 2d. What Table S2 does not contain — 🔴

| Item | Status |
|---|---|
| bead attachment ("5'-Bead-Linker-") | bead material and linker chemistry not given 🔴 |
| Tn5 adaptor | not listed 🔴; inferred in §3 step 8 from the printed N7xx primer and a later same-lab statement 🟡 |
| Index 1 and Read 2 sequencing primers | not listed 🔴; the stock Nextera primers are derived (§2c) 🟡 |
| i7 index orientation (synthesised bases vs index name) | not stated 🔴 |
| every reaction condition | in the unfetched STAR Methods 🔴 |

An earlier revision of this note carried an item-by-item comparison against the
scg_lib_structs drawn page and treated the agreements as confirmation. Agreement with that
page is not evidence, so the comparison is gone; what remains is what Tables S1–S2 print
and what is derived from them.

## 3. Step by step

Order 🟡 *(derived from the oligo architecture of Tables S1–S2: which oligo can only act
after which, see §2c)*; quantities, temperatures and cycle numbers 🔴 (STAR Methods not
obtained).

Bead synthesis (before the experiment):

1. **Round 1**: A-1…A-96 attached to magnetic beads, one oligo per well of a 96-well plate
   (attachment chemistry 🔴). Bead strand:
   `bead · head · SMART_HANDLE · ACGT · bc1 · L1`.
2. **Pool, split into B-1…B-96**, primer extension: B anneals by `L1'`; the bead strand is
   extended to `… · L1 · bc2 · L2`.
3. **Pool, split into C-1…C-96**, extension: C anneals by `L2'`; bead strand becomes
   `… · L2 · bc3 · UMI · T30`.
4. **Strip the complementary strand**, so the beads carry single-stranded seqA oligos.
   🟡 *(required by the next step — oligo-dT capture needs a free single strand — but no
   source we hold gives the method or conditions 🔴.)*

Library:

5. **Load cells** into the agarose microwell array, then beads; **lyse**; mRNA poly(A)
   anneals to the bead T30. Beads are collected (magnet) and pooled. 🟡 *(the magnetic
   retrieval is 🟢 in the same lab's adrenal paper; the 2018 conditions are 🔴)*
6. **RT with template switching** on the pooled beads (MMLV-type RT + TSO LNA): first
   strand `bead · … · UMI · T30 · cDNA · CCC`, then the TSO's `GGG` pairs with the
   untemplated `CCC` and the RT copies the TSO, adding `GAAT · SMART_HANDLE'` to the
   3' end. 🟡 ([template switching](../ref/concepts/template-switching.md))
7. **Whole-transcriptome PCR with TSO-PCR alone.** Both ends carry the SMART handle, so one
   primer amplifies full-length cDNA; the bead-side 20-nt head (with the I-SceI site) is
   copied at most once (first second-strand synthesis runs to the oligo's 5' end), but the
   primer anneals at the handle, so exponentially amplified products start at the handle. A
   single primer on two identical ends is nanoCAGE's **semi-suppressive PCR**: when both ends carry the same sequence, short
   products fold into intramolecular panhandles and are amplified poorly, while long cDNA
   amplifies (nanoCAGE Fig. 1b, read 🟢 for the principle). Purify.
8. **Tagment** the purified cDNA. The Guo-lab web protocol (§9) says "Nextera XT
   tagmentation reactions (Illumina, cat#FC-131-1096)" with the kit's N5xx primers
   replaced by the P5 primer 🟢 (lab protocol, undated). The same lab's 2020 adrenal paper
   instead names: a "customized
   transposase" from the **Vazyme TruePrep DNA Library Prep Kit V2 (TD513)**, "which
   carries two identical insertion sequences" (🟢 adrenal paper, 2020 — for 2018 🟡).
   Under either reading, the only Tn5-side primer is an N7xx (s7) primer, so the library's
   cDNA-side end is **s7-ME** (`nextera.ADAPTOR_S7`). With Nextera XT the s5-ME
   insertions present in that kit give no product with these primers. Which enzyme the
   2018 runs used is 🔴 (STAR Methods unread).
   🟡 (inferred; adapter bases not printed in any source held 🔴)
9. **Library PCR with P5 + one N7xx.** Three fragment classes 🟡 *(computed from the two
   printed primers and the inferred symmetric adaptor)*:
   - **bead end** (`SMART_HANDLE · AC · barcodes · UMI · T30 · insert · ME' · s7'`):
     P5 primes the handle end, N7xx the s7 end → **the library**;
   - **internal** (`s7 · ME · insert · ME' · s7'`): same adapter both ends, panhandle
     suppression and no P5 site → lost;
   - **TSO end** (`s7 · ME · insert · CCC · ATTC · SMART_HANDLE'`): has a handle, but P5's
     3' `AC` mismatches the `GA` that follows the handle there (§2c) → not amplified, or
     poorly.
10. Size-select / sequence (🔴 conditions).

## 4. Final library — 🟡 *(assembled from the Table S1/S2 oligos)*

```
5'- P5 · GCCTGTCCGCGG · SMART_HANDLE · ACGT · bc1(6) · L1(15) · bc2(6) · L2(15) · bc3(6) · UMI(6) · T30 · <cDNA, 3' end of transcript> · ME' · s7' · i7'(8) · P7' -3'
```

- 5' adapter up to the poly(T): 122 nt; 3' adapter (`ME' · s7' · i7' · P7'`) 66 nt, both
  assembled from the Table S2 primers and `lib/` constants 🟡 *(computed)*. The `ME' · s7'`
  part depends on the inferred symmetric adaptor (§3 step 8) and is drawn dotted.
- seqB beads: drop the `GT`; everything else the same.
- **No i5 index**: the P5 side carries the 12-nt spacer and handle, not an index.
  Libraries are single-indexed by i7.

## 5. Read layout / sequencing

| Read | Primer | Reads | Length needed |
|---|---|---|---|
| Read 1 | custom `Read1 SeqA-GT` (or `SeqB-TAC`) 🟢 | bc1 · L1 · bc2 · L2 · bc3 · UMI, then poly(T) | ≥ 54 nt to reach the end of the UMI 🟡 *(computed)* |
| Index 1 | `nextera.INDEX1_PRIMER` 🟡 *(derived, §2c)* | i7, 8 nt (the reverse complement of the index as Table S2 writes it) | 8 |
| Read 2 | `nextera.READ2_PRIMER` 🟡 *(derived, §2c)* | cDNA, from the Tn5 cut towards the poly(A) | — |

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
- 🟡 Tn5 adapter at the library end = s7-ME: follows from the printed N7xx primers. The
  undated lab web protocol names Nextera XT, while the 2020 adrenal paper names a TruePrep
  enzyme with "two identical insertion sequences". Neither is stated for 2018. This is
  the one inferred region on the public page, and it carries all three required signals
  (page caveat, `INFERRED — …` step caption, dotted bases).
- 🟡 `L1 = CGACTCACTACAGGG` differs from the T7 promoter tail `CGACTCACTATAGGG` by one
  base — probably a coincidence of design; not used as a promoter here.

## 7. How this note was made (tool evaluation)

`tools/get_sources.py` found no open copy of the paper (Cell/ScienceDirect return a
Cloudflare challenge). The two supplementary tables were first fetched by hand from the
scg_lib_structs repository's `data/Microwell-seq/` mirror of Elsevier's `mmc1`/`mmc2`;
`catalogue/extra_sources.tsv` now declares that mirror with the `scg:` handler, so the
tool fetches them itself and they land as `scg_*` where
`tools/selftest_source_oligos.py` can see them. Three same-lab papers plus the nanoCAGE
paper came from Europe PMC. `tools/doctext.py` copied the Europe PMC XML through as raw XML (tags kept);
the `.xml.txt` twins were re-made by stripping tags. `tools/scrape_primers.py` found all
Table S2 and Table S1 rows, but joined the one-sequence-per-line barcode text files into a
single 5472-nt "sequence".

## 8. Bounded pass 2026-10-17 (round breadth-first-12)

- Defining *Cell* article re-tried: `get_sources.py` still records it as a manual row
  (no open full text); a direct cell.com PDF request returned a Cloudflare 403 challenge,
  and Europe PMC core metadata for PMID 29474909 reports `inPMC=N`, `isOpenAccess=N`
  (full-text links: cell.com PDF and the DOI only). STAR Methods remain unread 🔴.
- Supplementary Table S2 (`upstream_1-s2.0-S0092867418301168-mmc2.xlsx.txt`, lines 1–9
  and the N701–N819 rows) re-compared with the model: bead seqA/seqB, TSO LNA, TSO-PCR,
  P5 `…GCAGAGT-s-A-s-C`, Read1 SeqA-GT/SeqB-TAC and N701 agree with the declared
  segments 🟢. No sequence change.
- Model correction: the graph's sequencing plan previously declared acquired cycle
  counts (Read 1 54, Index 1 8). Those are informative minima computed from the oligos
  (§5), not published run lengths, so the plan now declares none and states why; the page
  sequencing intro says the same.
- The same-lab adrenal paper (`same-lab_PMC7396412.xml.txt`, lines 3–24) only reports
  BGI MGISEQ2000 versus Illumina HiSeq sequencing of a later library; it does not state
  2018 read lengths.

## 9. Follow-up pass 2026-10-17 (round deferred-followup-02)

- Defining *Cell* article re-tried by every legitimate route. Unpaywall/OpenAlex: bronze
  open access at the publisher only (`https://www.cell.com/article/S0092867418301168/pdf`),
  no repository copy. cell.com fulltext/showPdf and ScienceDirect return Cloudflare 403
  challenge pages. Not in PMC or Europe PMC. Manual URL:
  http://www.cell.com/article/S0092867418301168/pdf. STAR Methods remain unread 🔴.
- GEO GSE108097 series and sample GSM2889480 (text view): the extraction protocol says only
  "Library preparations were performed with the Microwell-seq protocol". Instrument: HiSeq
  2500. SRA SRX3471805 run SRR6377012 has 2 × 150 nt deposited spots. These are deposit
  records, not a stated read design.
- **New original found: the Guo lab's own protocol** on the Mouse Cell Atlas website
  (https://bis.zju.edu.cn/MCA/, section "Microwell-seq protocol";
  `lab_MCA_microwell-seq_protocol_page.html.txt` lines 175–252) plus the linked docx
  `data/Protocol_for_making_microwells_and_barcoded_beads.docx`. Neither is dated or tied
  to a version, so applicability to the 2018 runs is 🟡. Read 🟢 (lab protocol):
  - beads: 28 µm × 35 µm microwells; carboxyl magnetic beads (20–25 µm); round-1
    amino-oligos coupled with EDC; rounds 2 and 3 by Phanta Master Mix (Vazyme) extension
    (5 cycles, then 1 cycle); complementary strands stripped at 95 °C (docx).
  - lysis 0.1 M Tris pH 7.5, 0.5 M LiCl, 1% LiSDS, 10 mM EDTA, 5 mM DTT, 10 min (line 193).
  - RT: SuperScript II, 1 M betaine, 6 mM MgCl2, 1 µM TSO_RT primer, 42 °C 90 min (line 199).
  - Exonuclease I (NEB), 1 U/µl, 37 °C 45 min, "to remove oligonucleotides which did not
    capture mRNA" (line 205).
  - WTA: KAPA HiFi HotStart ReadyMix + 0.4 µM TSO_PCR, 14 cycles (line 209).
  - "Nextera XT tagmentation … We replace the index 2 primers (N5xx) in the kit with our P5
    primer", 13–15 cycles, 0.6× AMPure, 400–700 bp peak; "sequenced on the Illumina
    NextSeq 500 using … Read1CustSeq" (line 213).
  - Table 1 (lines 229–252): TSO_RT `AAGCAGTGGTATCAACGCAGAGTGAAT/rG/rG/iXNA_G`, TSO_PCR,
    P5 `…GCAGAGT*A*C` and Read1CustSeq `GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTACGT`. All
    four agree with the Table S2-derived model (checked against `microwell_seq.TSO`,
    `P5_PRIMER` and `CUSTOM_R1`). Scraper run on both files: those 4 sequences only.
- Model changes: the graph now names SuperScript II, Exo I as its own digestion step, the
  KAPA PCR, the Phanta extension and the Nextera XT enzyme, all with lab-protocol
  provenance. The Tn5 adapter is relabelled from "symmetric s7" to "s7-ME at the
  N7xx-selected end" (still inferred). The page caveat and sequencing intro now give the
  NextSeq 500 (lab page) versus HiSeq 2500 (GEO) instrument difference. Construct bases are
  unchanged.
- Still deferred: the *Cell* STAR Methods (whether the lab page matches the 2018 runs), the
  2018 read lengths, publisher mmc3–mmc7, and the bases of the Tn5 adapter.

## 10. STAR Methods reading (deferred follow-up 03, 2026-10-17)

Source: `j.cell.2018.02.001.pdf.txt` (Cell 172:1091–1107.e1–e7). Verified as the defining
article (printed DOI 10.1016/j.cell.2018.02.001, title, Han … Guo). Line numbers below
point into the `.txt` twin.

- 🟢 Key Resources Table (lines 2510–2565): bead seqA/seqB as in Table S2; "Primer for
  Reverse transcription: AAGCAGTGGTATCAACGCAGAGTGAAT/rG/rG/iXNA_G"; "Primer for cDNA
  amplification: AAGCAGTGGTATCAACGCAGAGT"; library/sequencing primers "see Table S2".
- 🟢 Bead synthesis (lines 2680–2700): 20–25 µm carboxyl magnetic beads; round 1
  EDC coupling of amino-modified A oligos (50 µM) in 0.1 M MES in a 96-well plate;
  rounds 2 and 3 by PCR in 1× Phanta Master Mix with 5 µM B or C oligos (round 2:
  94 °C 5 min, 5 cycles 94 °C 15 s / 48.8 °C 4 min / 72 °C 4 min; round 3: 48.8 °C 20 min,
  72 °C 4 min); complementary chains removed at 95 °C for 6 min with magnetic separation,
  twice. This confirms step 4's strand stripping (§3).
- 🟢 Cell collection and lysis (lines 2701–2709): 100,000 cells/ml, 1,000,000 beads/ml;
  doublets blown out by capillary; beads loaded on a magnet; cold lysis buffer (0.1 M
  Tris-HCl pH 7.5 per correction, 0.5 M LiCl, 1% SDS, 10 mM EDTA, 5 mM DTT) 12 min; 6× SSC
  and Tris washes; 50,000 beads per tube.
- 🟢 RT (lines 2710–2715): "the instructions from the Smart-seq2 protocol were followed";
  200 U SuperScript II, 1 M betaine, 6 mM MgCl2, 2.5 mM DTT, 1 mM dNTP, TSO; 42 °C 90 min;
  TE-SDS wash.
- 🟢 Exonuclease I (lines 2716–2720): 1 U/µl, 37 °C **60 min** "to remove oligonucleotides
  that did not capture mRNA" (the lab web page says 45 min).
- 🟢 cDNA amplification (lines 2721–2727): four tubes, KAPA HiFi HotStart ReadyMix,
  0.1 µM TSO_PCR; 98 °C 3 min; 4 cycles (65 °C anneal) + 10–14 cycles (67 °C anneal);
  AMPure XP.
- 🟢 Transposase fragmentation and selective PCR (lines 2729–2736): "fragmented by a
  customized transposase that carries two identical insertion sequences … from TruePrep
  DNA Library Prep Kit V2 for Illumina (Vazyme)"; "We replaced the index 2 primers (N5xx)
  in the kit with our P5 primer … to specifically amplify fragments that contain the 3′ end
  of transcripts. Other fragments will form self-loops, which impede their binding to PCR
  primers." PCR 72 °C 3 min, 98 °C 30 s, 14 cycles (60 °C), AMPure XP; 400–700 bp peak;
  "sequencing on an Illumina Hiseq system". So the 2018 enzyme is the TruePrep V2
  single-insertion-sequence transposase, not Nextera XT; the insertion-sequence bases are
  not printed and s7-ME stays inferred 🟡.
- 🟢 Data processing (lines 2816–2820): reads lacking both linkers or poly(T) filtered;
  cell barcode and UMI from Read One. No read lengths are printed.
- Model corrections: conditions and sources moved from the undated lab web page to the
  STAR Methods (Exo I 60 min, WTA 4 + 10–14 cycles, library PCR 14 cycles, lysis buffer,
  bead synthesis); tagmentation relabelled TruePrep V2 with identical insertion sequences;
  the former discarded "tagmentation waste" replaced by a PCR off-product, the self-looped
  internal fragments, left unresolved (not physically removed); page caveat and
  sequencing text updated. Earlier statements in §§2d, 3, 5 and 6 that conditions are
  "unfetched" are superseded by this section.
- Still open 🔴 (not confirmed absent): insertion-sequence bases of the Vazyme TruePrep V2
  transposase (vendor document not read); acquired read lengths; publisher supplements
  mmc3–mmc7 (not on disk; mmc1/mmc2 are the scg-mirrored publisher files); bead linker
  chemistry ("5′-Bead-Linker-"); i7 index orientation question in §2c.
