# scifi-RNA-seq — plate pre-indexing by RT, then droplet overloading with a ligated bead barcode

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. Relationships marked 🟡 *(computed)* were worked out with `lib/` while
> writing this note; they are not yet asserted in a self-test, because this protocol has
> no `tools/` module yet (status `notes` in `catalogue/ours.tsv`). Claims taken only from
> the upstream scg_lib_structs page are 🟡 *(upstream)*: a secondary source.

**scifi-RNA-seq** (single-cell combinatorial fluidic indexing RNA-seq) —

- Preprint (defining, used for the directory name): Datlinger P\*, Rendeiro AF\*,
  Boenke T, Krausgruber T, Barreca D, Bock C. "Ultra-high throughput single-cell RNA
  sequencing by combinatorial fluidic indexing." *bioRxiv* 2019.
  doi:[10.1101/2019.12.17.879304](https://doi.org/10.1101/2019.12.17.879304) (v1, posted
  18 Dec 2019; "all rights reserved, no reuse allowed").
- Journal version: "Ultra-high-throughput single-cell RNA sequencing and perturbation
  screening with combinatorial fluidic indexing." *Nature Methods* 2021.
  doi:[10.1038/s41592-021-01153-z](https://doi.org/10.1038/s41592-021-01153-z),
  PMID 34059827, PMC7612019. Adds a single-cell CRISPR screen (T-cell receptor
  activation) and the step-by-step protocol (v2021-01-25) as a supplement.
- Code: github.com/epigen/scifiRNAseq (preprint, Data availability).
- Related method in the catalogue: **scifi-ATAC-seq** (doi:10.1101/2023.09.17.558155;
  *Genome Biology* 2024, doi:10.1186/s13059-024-03235-5) — same scifi idea on chromatin;
  separate directory, not covered here.

Sources read (fetched by `tools/get_sources.py`, into
`$CHEM_DATA/sources/scifi-rna-seq__10.1101+2019.12.17.879304/`, never committed):

| File | What | Used for |
|---|---|---|
| `2019.12.17.879304_v1.full.pdf` | the preprint, incl. Online Methods and Supplementary Figures | principle (main text), **all reaction conditions** (Online Methods), Supp. Fig. 2 (sequence design — an image, read visually; see below) |
| `2019.12.17.879304_media-1.xlsx` | preprint Supplementary Table 1 | **every oligo** (sheet "scifi-RNA-seq oligos"), the 384 RT primers (sheet "384 RT primers") |
| `2019.12.17.879304_media-2.xlsx` | sequencing performance per sample | not chemistry |
| `s41592-021-01153-z_..._MOESM2_ESM.xlsx` | journal Supplementary Table (oligos) | same as media-1 plus **P7-009 … P7-016** |
| `s41592-021-01153-z_..._MOESM5_ESM.pdf` | journal supplement: **scifi-RNA-seq step-by-step protocol v2021-01-25** | master-mix tables, scATAC v1.0 *and* v1.1 (Next GEM) variants, Nextera XT alternative, NovaSeq read structure, Tn5 assembly |
| `s41592-021-01153-z_..._MOESM6_ESM.xlsx` | CRISPR-screen table | checked that its 13-nt well barcodes are the RT-primer barcodes |
| `s41592-021-01153-z_..._MOESM3/4/7/8_ESM.xlsx` | sequencing metrics, datasets, nuclei/bead counts, droplet sizes | not chemistry |
| `upstream_scifi-RNA-seq.html` | scg_lib_structs page (a redrawing of Supp. Fig. 2) | second source, checked below |
| `sfig2/p20-000.png` (+ crops) | Supp. Fig. 2 extracted from the preprint PDF page 20 with `pdfimages` | the drawn construct at every step, both sequencer layouts |

Not obtained:

| What | Why | URL |
|---|---|---|
| Journal full text (main text + Methods of the Nature Methods version) | PMC returned a reCAPTCHA page, which was saved as `s41592-021-01153-z_PMC7612019.html` | https://pmc.ncbi.nlm.nih.gov/articles/PMC7612019/ · https://www.nature.com/articles/s41592-021-01153-z |
| `MOESM1_ESM.pdf` contents | image-only PDF (title "Bock RS.pdf", printed to PDF — evidently the Reporting Summary); `.txt` twin is empty | (on disk; OCR would be needed, likely not chemistry) |

---

## 1. What it is

Two rounds of cell barcoding, the second one done in a commercial droplet instrument that
is deliberately **overloaded** 🟢 (preprint main text, Online Methods):

1. **Round 1 — plate.** Fixed/permeabilised cells or nuclei are distributed into a
   384-well plate (10,000 per well) and reverse transcribed *in situ* with a
   **well-specific barcoded oligo-dT primer** that also carries an 8-nt UMI and the
   TruSeq Read 1 site, and a **5′ phosphate**. Maxima H Minus adds the untemplated C's at
   the end of the cDNA.
2. **Round 2 — droplets.** All wells are pooled and up to 765,000 pre-indexed nuclei
   (1.53 million tested for Next GEM) are loaded into **one** Chromium channel with
   **10x scATAC gel beads** (not the scRNA beads). Most droplets hold several cells; the
   cell is resolved by the *pair* (round-1 well barcode, round-2 bead barcode). The bead
   oligo is joined to the RT primer's 5′ end by **ligation on a 3′-blocked bridge oligo**
   with Ampligase, cycled 12× between denaturation and 59 °C.
3. **Bulk.** Emulsion broken, cDNA purified, then **template switching in bulk** on the
   pooled single-stranded cDNA, cDNA PCR, tagmentation with a custom **i7-only Tn5**, and
   PCR with partial P5 + an indexed P7 primer.

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | **Droplet overloading** made usable by a plate pre-index (combinatorial fluidic indexing) | sci-RNA-seq does all rounds in plates; 10x scRNA-seq has one cell per droplet |
| 2 | Bead barcode attached by **templated thermoligation inside droplets** (bridge oligo + Ampligase, 5′-phosphate on the RT primer) instead of by priming from the bead | split-pool methods ligate barcodes in wells (SPLiT-seq, sci-RNA-seq3) |
| 3 | Uses the **scATAC** bead (P5 · bc16 · s5), so the bead barcode is read as **Index 2** | 10x scATAC-seq itself; scifi-ATAC-seq |
| 4 | **Template switching after** the emulsion, in bulk, on purified first-strand cDNA (the C's were added in situ) | normally TSO is present during RT — [template switching](../ref/concepts/template-switching.md) |
| 5 | **i7-only transposome** (Tn5 loaded with s7-ME only), so every Tn5 end is P7-side and only bead-end fragments amplify | [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md); Drop-seq/10x 3′ use Nextera with a P5 bead primer for the same effect |

Builds on: [reverse transcription](../ref/concepts/reverse-transcription.md) (barcoded
anchored oligo-dT), [template switching](../ref/concepts/template-switching.md) (SMART
handle TSO, `rGrGrG`), [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md).

## 2. Oligos

🟢 Verbatim from Supplementary Table 1 (`media-1.xlsx` / journal `MOESM2_ESM.xlsx`, sheet
"scifi-RNA-seq oligos"), modifications as written there. The table's legend: `[5Phos]`
5′ phosphorylation, `[ddC]` dideoxy-C "to block 3′ extension", `r` = next base is RNA,
`N` any base, `V` not T. (The legend also says the mosaic end is "underlined"; the
formatting is lost in the text twin.) The preprint and journal tables are **identical**
for every oligo they share 🟡 (computed: all 16 named oligos and all 384 RT primers
compared string by string); the journal table adds P7-009 … P7-016.

```
Tn5-top_ME                     [5Phos]CTGTCTCTTATACACATCT
Tn5-bottom_Read2N              GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG
pUC19-FWD                      AAGTGCCACCTGACGTCTAAG              (Tn5 activity test only)
pUC19-REV                      CAACAATTAATAGACTGGATGGAGGCGG       (Tn5 activity test only)
Bridge-Oligo_truseq_ddC        CGTCGTGTAGGGAAAGAGTGTGACGCTGCCGACGA[ddC]
Template-Switching-Oligo(TSO)  AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG
Partial-P5                     AATGATACGGCGACCACCGAGA
TSO_enrichment_primer          AAGCAGTGGTATCAACGCAGAGT
```

**Round-1 RT primers**, sheet "384 RT primers" (4 plates × 96 = 384, `SCIFI_LIG384_001`
… `_384`, index IDs `SASI_001` … `SASI_384`). First three as written:

```
SCIFI_LIG384_001  [5'phos]ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNNNAAGTGATTAGCAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN   index AGTGATTAGCA  fixed A / A
SCIFI_LIG384_002  [5'phos]ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNNNAGAATCCCCCTAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN   index GAATCCCCCTA  fixed A / A
SCIFI_LIG384_003  [5'phos]ACACTCTTTCCCTACACGACGCTCTTCCGATCTNNNNNNNNACCTGGGAAACTATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN   index CCTGGGAAACT  fixed A / A
```

(This sheet writes the phosphate `[5'phos]`, the oligo sheet `[5Phos]` — same thing.)

**Indexed P7 primers** (column `P7_index_seq` = the Illumina name's i7 bases):

```
P7-001  CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGG  NexteraXTv2_N701  TAAGGCGA
P7-002  CAAGCAGAAGACGGCATACGAGATCTAGTACGGTCTCGTGGGCTCGG  NexteraXTv2_N702  CGTACTAG
P7-003  CAAGCAGAAGACGGCATACGAGATTTCTGCCTGTCTCGTGGGCTCGG  NexteraXTv2_N703  AGGCAGAA
P7-004  CAAGCAGAAGACGGCATACGAGATGCTCAGGAGTCTCGTGGGCTCGG  NexteraXTv2_N704  TCCTGAGC
P7-005  CAAGCAGAAGACGGCATACGAGATAGGAGTCCGTCTCGTGGGCTCGG  NexteraXTv2_N705  GGACTCCT
P7-006  CAAGCAGAAGACGGCATACGAGATCATGCCTAGTCTCGTGGGCTCGG  NexteraXTv2_N706  TAGGCATG
P7-007  CAAGCAGAAGACGGCATACGAGATGTAGAGAGGTCTCGTGGGCTCGG  NexteraXTv2_N707  CTCTCTAC
P7-008  CAAGCAGAAGACGGCATACGAGATCAGCCTCGGTCTCGTGGGCTCGG  NexteraXTv2_N710  CGAGGCTG
P7-009  CAAGCAGAAGACGGCATACGAGATAAGAGGCAGTCTCGTGGGCTCGG  NexteraXTv2_N711  AAGAGGCA   (journal table only)
P7-010  CAAGCAGAAGACGGCATACGAGATGTAGAGGAGTCTCGTGGGCTCGG  NexteraXTv2_N712  GTAGAGGA   (journal table only)
P7-011  CAAGCAGAAGACGGCATACGAGATGCTCATGAGTCTCGTGGGCTCGG  NexteraXTv2_N714  GCTCATGA   (journal table only)
P7-012  CAAGCAGAAGACGGCATACGAGATATCTCAGGGTCTCGTGGGCTCGG  NexteraXTv2_N715  ATCTCAGG   (journal table only)
P7-013  CAAGCAGAAGACGGCATACGAGATACTCGCTAGTCTCGTGGGCTCGG  NexteraXTv2_N716  ACTCGCTA   (journal table only)
P7-014  CAAGCAGAAGACGGCATACGAGATGGAGCTACGTCTCGTGGGCTCGG  NexteraXTv2_N718  GGAGCTAC   (journal table only)
P7-015  CAAGCAGAAGACGGCATACGAGATGCGTAGTAGTCTCGTGGGCTCGG  NexteraXTv2_N719  GCGTAGTA   (journal table only)
P7-016  CAAGCAGAAGACGGCATACGAGATCGGAGCCTGTCTCGTGGGCTCGG  NexteraXTv2_N720  CGGAGCCT   (journal table only)
```

P7-001 … P7-008 are also printed in the step-by-step protocol (MOESM5, step 10A) with
the same sequences 🟢.

**Not in any table — the 10x scATAC gel-bead oligo.** Shown only in Supp. Fig. 2
(an image; "Gel Bead Oligo Primer with 737k possible microfluidic barcodes (16b)") and on
the upstream page, which agree 🟡 (read from the figure image + upstream; not
text-verifiable):

```
gel bead oligo   bead-|-5'- AATGATACGGCGACCACCGAGATCTACAC <16-nt bead barcode> TCGTCGGCAGCGTC -3'
```

### How the oligos interlock — 🟡 (computed)

**RT primer** (86 nt with N's, all 384 checked):

```
5'-phos · illumina.TRUSEQ_READ1 (33) · UMI N8 · fixed1 (1) · SASI index (11) · fixed2 (1) · T30 · V · N -3'
```

- Every one of the 384 sequences equals exactly that concatenation, with the table's
  own `index_seq`, `fixed_base1_N` and `fixed_base2_V` columns. The 11-nt indices are
  all distinct, minimum pairwise Hamming distance **5**.
- The two **fixed bases** follow the same layout on all four plates. fixed1 is set
  per plate row: A (rows A–B), T (C–D), C (E–F), G (G–H) — 96 each. fixed2 is *not*
  per row: it runs in blocks of 32 wells in row-major order — A for A01–C08, C for
  C09–F04, G for F05–H12 — 128 each. Supp. Fig. 2 labels them "Fixed bases", Read 1 as "UMI (8) + fixed N (1) +
  Round1 (11) + fixed V (1)" (figure text). Fixed2 being V presumably stops the
  barcode from running into the T30 🟡; why fixed1 is varied by row is not stated 🔴.
- The 13-nt "scifi-RNA-seq barcodes" in the CRISPR-screen table (MOESM6) are
  fixed1 + index + fixed2: all 372 listed occur among the 384 primers.

**Bridge oligo ↔ bead oligo ↔ RT primer.** revcomp(bridge without ddC) =
`TCGTCGGCAGCGTC` + `ACACTCTTTCCCTACACGACG`, i.e.

- its 3′ part (14 nt) pairs with the **whole s5** (`nextera.S5`) at the 3′ end of the
  gel-bead oligo, and
- its 5′ part (21 nt) pairs with the **first 21 nt of TruSeq Read 1**
  (`illumina.TRUSEQ_READ1[:21]`) at the 5′ end of the RT primer,

so the bead's 3′-OH and the RT primer's 5′-phosphate abut at a nick on the bridge:
a 14 + 21 split. The 3′ `ddC` lies one base beyond s5, opposite the last base of the
16-nt bead barcode (variable, so it pairs only by chance); its stated job is to stop the
bridge from being extended 🟢 (legend). Supp. Fig. 2 draws the bridge in exactly this
register.

After ligation the junction reads `… <bead bc16> · TCGTCGGCAGCGTC · ACACTCTTTCCC…` —
s5 directly followed by TruSeq Read 1, no gap.

**TSO and its primer.** TSO = **`rt.SMART_HANDLE`** (23 nt) + `GAAT` + `rGrGrG` (30 nt;
all-DNA except the three ribo-G's, no LNA). The TSO_enrichment_primer **is**
`rt.SMART_HANDLE`. The `GAAT` spacer is not explained 🔴. The cDNA 3′ end after
switching is revcomp(TSO) = `CCCATTCACTCTGCGTTGATACCACTGCTT` (matches Supp. Fig. 2).

**Partial-P5** = `illumina.P5[:22]` — it stops 7 nt short of the bead barcode, so it
"retains round1 and round2 cell barcodes" (figure) and does not need to cover the bead
barcode.

**Tn5 oligos.** Tn5-top_ME = `nextera.ME_RC` (the 5′-phosphorylated, transferred-strand
partner), Tn5-bottom_Read2N = `nextera.ADAPTOR_S7` (s7 + ME). There is **no s5
adaptor** — the transposome is i7-only.

**P7 primers** = `illumina.P7` + 8 nt + `nextera.S7` (all 47 nt; the form built by
`nextera.n7xx_primer`, which inserts its argument as written — so P7-001 … 008 equal
`n7xx_primer(revcomp(P7_index_seq))` and P7-009 … 016 equal `n7xx_primer(P7_index_seq)`).
The s7 tail is the 15-nt s7 only, no ME bases. Index orientation is **inconsistent in
the table**:

- **P7-001 … P7-008**: the primer contains the **reverse complement** of the listed
  `P7_index_seq` (e.g. P7-001 contains `TCGCCTTA` for N701 `TAAGGCGA`) — the usual
  Illumina N7xx convention, so Index 1 reads the listed bases.
- **P7-009 … P7-016**: the primer contains the listed bases **as written** (P7-009
  contains `AAGAGGCA` for N711 `AAGAGGCA`). If synthesised as printed, Index 1 would read
  the reverse complement of the listed i7 for these eight. Either a transcription error
  in the journal table or a real, unflagged difference — see §7.

## 3. Step by step

Volumes and conditions 🟢 from the preprint Online Methods and the step-by-step protocol
(MOESM5); where they differ, both are given.

1. **Cells / nuclei.** Four entry routes in the protocol: methanol-fixed whole cells
   (1A), fresh nuclei (1B: NP-40/IGEPAL, digitonin, Tween-20), formaldehyde-fixed and
   snap-frozen nuclei, re-permeabilised after thawing (1C, "for primary cells"). Dilute
   to 5,000/µL.
2. **Round-1 RT, in situ** (384-well plate, 1 µL of 25 µM barcoded RT primer pre-dispensed
   per well). 10,000 cells/nuclei (2 µL) per well, 55 °C 5 min → ice (RNA secondary
   structure), + 7 µL mix: 1× RT buffer, 5 mM DTT, 0.5 mM dNTPs, 20 U RNaseOUT, 100 U
   **Maxima H Minus**. 50 °C 10 min; 3 × (8 °C 12 s, 15 °C 45 s, 20 °C 45 s, 30 °C 30 s,
   42 °C 2 min, 50 °C 3 min); 50 °C 5 min. **No TSO here.** Product, inside each nucleus:
   `5'-phos-R1-UMI-bc-T30VN-…first-strand cDNA…-CCC-3'` on its mRNA; the CCC is the RT's
   untemplated addition (main text 🟢; drawn in Supp. Fig. 2); see
   [template switching](../ref/concepts/template-switching.md).
3. **Pool and count.** Wash wells with PBS–1 % BSA, pool, spin 500 rcf, resuspend in
   **1× Ampligase buffer**, filter, pellet to < 15 µL, count (Fuchs-Rosenthal).
4. **Round-2 thermoligation in droplets.** Up to 765,000 nuclei in 15 µL + thermoligation
   mix (Ampligase buffer, 230 U Ampligase, Reducing Agent B, **2.3 µL of 100 µM bridge
   oligo**) to 80 µL (v1.0, Chip E) or 75 µL (v1.1 Next GEM, Chip H), load with **scATAC
   gel beads**, run the Controller. Emulsion: **12 × (95 °C 30 s, 59 °C 2 min)** in the
   protocol — **98 °C 30 s** in the preprint Methods. Cells lyse in the droplet; the bead
   oligo is released (Reducing Agent B dissolves the bead 🟡) and is ligated to the RT
   primer's 5′-phosphate across the bridge 🟢 (mechanism, main text). Thermocycling
   denatures each product off the bridge so a bridge can template the next ligation 🟡.
   First-strand product:
   `5'-P5 · bc16 · s5 · TruSeq R1 · UMI · bc13 · T30VN · cDNA · CCC-3'`.
   Variant (Appendix H, v1.0 only): thermoligation in 10x Barcoding Reagent + NAD⁺ and
   MgCl₂ when the Ampligase-buffer emulsion is unstable.
5. **Break and clean.** Recovery Agent; Dynabeads MyOne Silane cleanup (or MinElute
   column, Appendix G, when genomic DNA clumps the beads); 1.0× SPRIselect, 20 µL.
6. **Template switching, bulk.** 20 µL sample + 1× RT buffer, 4 % Ficoll PM-400, dNTPs,
   RNase inhibitor, **1.25 µL of 100 µM TSO**, 500 U Maxima H Minus in 50 µL; 25 °C
   30 min, 42 °C 90 min. The TSO's rGrGrG pairs with the CCC and the cDNA is extended
   across the TSO, adding `ATTCACTCTGCGTTGATACCACTGCTT` after the CCC 🟡 (computed; drawn
   in Supp. Fig. 2). 1.0× AMPure.
7. **cDNA enrichment**: NEBNext HiFi 2×, 0.5 µM each **Partial-P5** and
   **TSO_enrichment_primer**, SYBR Green; 98 °C 30 s; cycle (98 °C 20 s, 65 °C 30 s,
   72 °C 3 min) until > 1000 RFU; 72 °C 5 min. The TSO primer makes the second strand;
   both then amplify 🟡. Cleanup 0.8× then 0.6× AMPure; Qubit; Bioanalyzer.
8. **Tagmentation**, 1 ng per reaction (6–15 reactions per sample):
   - **10A (default), i7-only Tn5**: 1× TAPS-MgCl₂ buffer, 10 % DMF, 1.25 µL transposome
     (1:4.5) in 25 µL; 55 °C 10 min; stop with 2.5 µL 1 % SDS; 1.0× AMPure.
   - **10B, Nextera XT** (Amplicon Tagment Mix, s5 *and* s7 adaptors): 55 °C 5 min,
     Neutralization Buffer. The protocol states this loses 50 % of fragments with
     comparable results; i7-only is "conceptually cleaner, cheaper and more versatile".
9. **Library PCR**: + barcoded P7 primer (5 µL of 10 µM in 10A; 2.5 µL in 10B) and
   Partial-P5 (0.5 µL of 100 µM in 10A), NEBNext HiFi (10A) or Nextera PCR mix (10B), SYBR; **72 °C 3 min
   (gap fill-in)**, 98 °C 30 s, cycle (98 °C 10 s, 65 °C 30 s, 72 °C 30 s) to > 4000 RFU;
   72 °C 5 min. Several P7 primers per sample are allowed for colour balance on 2-colour
   sequencers. 0.7× AMPure per well, pool per P7 index, 0.8× AMPure.

Transposome (Appendix A / Methods): anneal Tn5-top_ME + Tn5-bottom_Read2N, 45 µM each in
1× annealing buffer (95 °C 3 min, 70 °C 3 min, ramp 2 °C/min to 25 °C), dilute with
180 µL water; 20 µL + 20 µL glycerol + 10 µL **EZ-Tn5 (Lucigen TNP92110)**, 25 °C 30 min.
Activity is checked in a **negative** qPCR on a 1961-bp pUC19 amplicon, because
i7–i7 fragments are PCR-suppressed and cannot be measured positively 🟢.

### What survives tagmentation — 🟡 (computed / inferred)

With the i7-only Tn5, every Tn5 end carries `s7 · ME`. Of the pieces of one cDNA:

- **bead end + Tn5 end** (P5 · bc16 · s5 · R1 · UMI · bc13 · T30VN · cDNA · ME · s7) —
  Partial-P5 primes one end, the P7 primer the other: **the library**.
- **Tn5 + Tn5** (s7 on both ends) — only the P7 primer binds; suppression PCR, lost.
- **Tn5 + TSO end** — no primer for the SMART handle in this PCR, lost.

So Read 2 starts at a Tn5 site somewhere upstream of the poly(A) and the library is
3′-tag counting, as in 10x 3′ and Drop-seq. With Nextera XT (10B) half the bead-end
fragments get s5 instead of s7 at the Tn5 end and are lost — the protocol's "50 % loss".
The upstream page draws exactly these three products 🟡 (upstream), agreeing.

## 4. Final library — 🟡 (assembled from the oligos above)

```
5'- P5 (29) · <bead bc, 16> · s5 (14) · TruSeq R1 (33) · <UMI, 8> · <fixed N, 1> · <round-1 bc, 11> · <fixed V, 1> · T30 · V · N · <cDNA, sense of first strand = antisense of mRNA> · ME' (19, CTGTCTCTTATACACATCT) · s7' (15) · <i7, 8> · P7' (24) -3'
```

- Fixed (non-insert) length: 113 nt left of the T30, 32 nt `T30VN`, 66 nt right of the
  cDNA: **211 nt** plus insert (computed with `lib/`).
- The 16-nt bead barcode sits **in the i5 index position** (between P5 and s5), exactly
  where a dual-index Nextera library has its i5 — which is why it can be read as Index 2
  with standard primers on NovaSeq.

Upstream's final library agrees base for base in the fixed parts 🟡 (checked with
`grep`): it writes the 21 cycles of Read 1 as 20 N + `V`, i.e. it shows fixed2 explicitly
and lumps fixed1 into the barcode ("13-bp RT barcode"), and it drops the final `N` of
`VN` into the cDNA — both just notation.

## 5. Sequencing

🟢 (preprint Methods; protocol step 12): NovaSeq 6000, SP / S1 / S2 100-cycle kits,
"standard sequencing primers", **Read 1 21, Index 1 (i7) 8, Index 2 (i5) 16, Read 2
78**. Aim for ≥ 10,000 reads per cell; load 2.0 nM.

| Read | Primer | Reads | 🟡 (computed) |
|---|---|---|---|
| Read 1, 21 nt | TruSeq Read 1 (`illumina.TRUSEQ_READ1`) | UMI 8 + fixed 1 + round-1 bc 11 + fixed 1 | 8+1+11+1 = 21: Read 1 stops exactly before the T30 |
| Index 1, 8 nt | Nextera i7 primer (`nextera.INDEX1_PRIMER`) | sample index | |
| Index 2, 16 nt | standard i5 read (NovaSeq, forward-strand workflow, primed off P5) | **bead barcode (round 2)** | |
| Read 2, 78 nt | Nextera Read 2 (`nextera.READ2_PRIMER`) | cDNA from the Tn5 end | reads mRNA-sense sequence 🟡 |

**NextSeq 500 variant** (Supp. Fig. 2 only — a figure panel, image text 🟡): 75-cycle kit,
Read 2 **47** cycles, and Index 2 read with a **custom index-2 primer**, because the
NextSeq reads i5 off the other strand and there is no standard primer site
next to the bead barcode on that side. Upstream gives this primer as
`GGAAAGAGTGTGACGCTGCCGACGA` ("Bead barcode sequencing primer (i5)") 🟡 (upstream; drawn
3′→5′ in Supp. Fig. 2). Computed: it is the reverse complement of `s5` + the first 11 nt
of TruSeq Read 1, i.e. it anneals across the ligation junction and reads the bead
barcode as its reverse complement — and it is **bases 11–35 of the bridge oligo**, so
the bridge minus its first 10 nt and ddC doubles as this primer.

### Upstream page vs the paper — summary

| Item | Upstream | Paper / supplement | Verdict |
|---|---|---|---|
| RT primer | R1 · UMI8 · 13-bp RT barcode · T30VN, 5′ phos | same; 13 = fixed1 + 11-nt index + fixed2 | agree (upstream lumps the fixed bases) |
| Bridge, TSO, Partial-P5, TSO primer, Tn5 oligos | as Table 1 | Table 1 | agree |
| Gel bead oligo | P5 · bc16 · s5 | figure only | agree |
| Template switching after the emulsion | yes | yes | agree |
| Tagmentation | "Tn5 homodimer" with s7-ME, 3 products | i7-only transposome (default) or Nextera XT | agree for 10A; upstream omits the XT variant |
| Index primer | "Nextera N7 index primer" P7 · i7 · s7 | P7-001…016 | agree (upstream does not show the P7-009…016 orientation issue) |
| **Read 2 length** | **47 cycles** | **78** (NovaSeq, Methods + protocol); 47 only in the NextSeq panel | upstream shows the NextSeq layout |
| **Index 2 primer** | custom bead-barcode primer | **standard** primers on NovaSeq; custom only for NextSeq | same — upstream drew the NextSeq panel only |
| Thermoligation details, scATAC v1.1 support, CRISPR screen | not mentioned | protocol / journal | upstream predates or omits them |

## 6. Variants in the sources

- **scATAC v1.0 (Chip E) vs v1.1 Next GEM (Chip H)** — only volumes differ (loading 75 vs
  70 µL; mix 80 vs 75 µL; glycerol fills); chemistry identical from cleanup on 🟢.
- **i7-only Tn5 (10A) vs Nextera XT (10B)** — see §3.8.
- **NovaSeq vs NextSeq read-out** — see §5.
- **CRISPR screen** (journal): gRNAs are listed per round-1 well (four wells per gRNA ×
  stimulation); the source table maps wells to perturbations, i.e. the perturbation was
  encoded by the round-1 barcode 🟡 (inferred from the table's columns). How the
  perturbation was introduced and whether a gRNA-enrichment library exists is in the
  journal Methods, which were not obtained 🔴.

## 7. Open questions

- 🔴 **P7-009 … P7-016 orientation**: the journal table writes the listed i7 bases
  forward inside the primer, unlike P7-001 … 008. Error in the table or the actual
  synthesised oligos? Needs the journal Methods or the epigen/scifiRNAseq sample sheets.
- 🔴 Purpose of the **row-wise fixed1 base** and the **`GAAT`** spacer in the TSO.
- 🔴 Thermoligation denaturation **95 °C** (protocol v2021-01-25) vs **98 °C** (preprint).
- 🟡 That the mRNA is gone (or irrelevant) by the time of bulk template switching — the
  12 thermoligation cycles at ≥ 95 °C and the silane/SPRI cleanups suggest the cDNA is
  single-stranded; not stated.
- 🔴 Journal-version Methods not read (PMC captcha): any changes to oligos or conditions
  in the published version beyond the supplements are unchecked.

## 8. How this note was made (tool evaluation)

`tools/get_sources.py` fetched the preprint, both preprint supplements, all eight
journal supplements and the upstream page; the PMC full text came back as a reCAPTCHA
page (saved under the article's name, not flagged as `(manual)`), and the journal
Reporting Summary is an image-only PDF whose `.txt` twin is empty. `scrape_primers.py`
put Supplementary Table 1 first and recognised `nextera.ME_RC`, `ADAPTOR_S7` and
`rt.SMART_HANDLE`. Supp. Fig. 2 — the only source for the bead oligo, the bridge register
and the NextSeq layout — is a raster image inside the preprint PDF; it was extracted
with `pdfimages` into `sfig2/` and read visually.

`scrape_primers.py --find` located every 🟢 oligo except the three RT primers: it
reports 0 locations for any query containing `N`/`V`, although the strings occur
verbatim in both RT-primer sheets (checked with `grep`; their N-free 3′ halves are found).
