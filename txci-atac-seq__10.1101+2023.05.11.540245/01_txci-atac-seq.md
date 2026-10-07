# txci-ATAC-seq — Tn5 pre-indexing on a 96-well plate, then overloading the 10x scATAC chip

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. Relationships marked 🟡 *(computed)* were worked out with `lib/` while
> writing this note; they are not yet asserted in a self-test, because this protocol has
> no `tools/` module yet (status `notes` in `catalogue/ours.tsv`).

**txci-ATAC-seq** (10X-compatible combinatorial indexing ATAC-seq) — Zhang H\*, Mulqueen RM\*,
Iannuzo N, Farrera DO, Polverino F, Galligan JJ, Ledford JG, Adey AC, Cusanovich DA.

- Preprint: "txci-ATAC-seq, a massive-scale single-cell technique to profile chromatin
  accessibility." *bioRxiv* 2023. doi:[10.1101/2023.05.11.540245](https://doi.org/10.1101/2023.05.11.540245)
  (v1, posted 14 May 2023; CC-BY-NC-ND 4.0).
- Published: same title (with a colon), *Genome Biology* 25:78 (22 Mar 2024).
  doi:[10.1186/s13059-023-03150-1](https://doi.org/10.1186/s13059-023-03150-1), PMC10958877.
  Not in the catalogue row; found by hand (Crossref title match) 🟡.
- Builds on: **dsciATAC-seq** (Lareau et al. 2019, Bio-Rad ddSEQ, Tn5-indexed then
  droplet-indexed) and the authors' own **sci-ATAC-seq**; the droplet part is the
  unmodified **10x Chromium Single Cell ATAC v1.1** kit (user guide CG000209 Rev D).
- The catalogue's "cited" row, nanoCAGE / semi-suppressive PCR (Plessy et al. 2010,
  doi:10.1038/nmeth.1470), is cited by upstream only to explain why A/A and B/B
  tagmentation products do not amplify; it is not a txci-ATAC paper.
- Protocol: protocols.io "TEN (10X)-compatible combinatorial indexing ATAC-seq"
  (`https://www.protocols.io/view/ten-10-x-compatible-combinatorial-indexing-atac-se-c3urynv6.html`)
  — **not fetched** (🔴 here).

Sources read (in `$CHEM_DATA/sources/txci-atac-seq__10.1101+2023.05.11.540245/`, never committed):

| File | What | Used for |
|---|---|---|
| `2023.05.11.540245_v1.full.pdf.txt` | preprint, incl. Supplementary Tables 4, 7, 8 printed at the end (lines ~2050–2215) | methods, oligos, cycle numbers |
| `PMC10958877.xml.txt` | Genome Biology full text (Europe PMC) | checking what changed on publication |
| `13059_2023_3150_MOESM5_ESM.xlsx.txt` | Additional file 5: Tables S4–S8 | **every oligo**, the 96 Tn5 barcodes, plate design |
| `13059_2023_3150_MOESM1_ESM.docx.txt` | Additional file 1: supplementary figure legends | Fig. S1/S2/S18 (barcode swapping, PCR efficiency) |
| `13059_2023_3150_MOESM6_ESM.docx.txt` | Additional file 6: review history | skimmed; nothing new on chemistry |
| `upstream_txci-ATAC-seq.html.txt` | scg_lib_structs page | second source, checked below |
| `upstream_Tn5ME-B_barcode.xlsx.txt` | upstream's copy of the Tn5 barcode plate | compared with Table S5 |

Not fetched: the protocols.io page (above). Additional files 2–4 (differential peaks,
motifs, KEGG) were not needed — no chemistry.

---

## 1. What it is

Two-level combinatorial indexing for scATAC on a stock 10x Chromium ATAC chip 🟢:

1. **Pre-index** nuclei in bulk: each well of a 96-well plate holds a Tn5 loaded with a
   **different 8-nt barcode** (Illumina "iTSM" plate, gift of Illumina). Each sample gets
   its own set of wells — this is also the **sample multiplexing**.
2. **Pool and overload** the 10x chip: 75,000–200,000 nuclei in one lane instead of the
   recommended maximum of 15,300. Many droplets hold several nuclei; they are told apart by
   their Tn5 barcode.
3. In-droplet PCR adds the 16-nt **GEM (bead) barcode**; Sample Index PCR adds an 8-nt
   **i7** per 10x lane.

Cell = i7 × GEM barcode × Tn5 barcode (concatenated in that order in the analysis) 🟢.
Up to ~22-fold more nuclei per lane than the standard workflow at the same collision rate 🟢.

| | New thing here | Concept / where else |
|---|---|---|
| 1 | Tn5 adaptor **B carries a TruSeq Read 2 handle + 8-nt barcode** instead of the Nextera s7, so the stock 10x bead primer (P5-GEM-s5) still primes the A end and **standard Illumina read primers** read everything — no custom sequencing recipe (the contrast with dsciATAC's custom primers is ours 🟡) | [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md) |
| 2 | **Short SBS primer spiked into the GEM** turns the 10x in-droplet *linear* amplification into *exponential* PCR, which suppressed in-droplet Tn5-barcode swapping (collision 46.0 % → 6.6 % in the cell-line true barnyard) | blocking oligo and decoy DNA were the two losing alternatives |
| 3 | **10 % / 90 % GEM split**: a tenth of the emulsion is processed first as a QC library | — |
| 4 | **Phased-txci-ATAC-seq** (preprint name: **Fast-txci-ATAC-seq**): nuclei frozen in tube strips per sample, tagmented straight after thawing | §4 |

The barcode swapping it fixes 🟢 (main text, Fig. S1c): in an overloaded droplet, restricting
to GEMs with a single Tn5 barcode dropped the collision rate from ~40–47 % to 3–9 %, so
most "doublets" were not two nuclei but Tn5 barcodes moving between fragments in the
droplet. The mechanism is drawn in Fig. S1c (not read; figure only) 🔴; the paper's stated
rationale for the blocking oligo is that free Tn5 adaptors act as primers 🟢.

## 2. Oligos

🟢 Verbatim from Additional file 5 (Tables S4, S7, S8; identical in the preprint's printed
tables). Modification notation as the source writes it: `[phos]` = 5' phosphate,
`/ideoxyU/` = internal 2'-deoxyuridine, `/3InvdT/` = 3' inverted dT (the main text calls it an inverted dideoxythymidine); lowercase in Full
SBS = the barcode (source's own note).

Tn5 linker oligos (Table S4):

```
Tn5ME-A    TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
Tn5ME-B    CGTGTGCTCTTCCGATCTNNNNNNNNAGATGTGTATAAGAGACAG      N = 8-nt Tn5 barcode (Table S5)
Tn5MErev   [phos]CTGTCTCTTATACACATCT
```

TruSeq i7 index primers for Sample Index PCR (Table S7), with the "Barcode" column:

```
P7.S701  CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT   TAAGGCGA
P7.S702  CAAGCAGAAGACGGCATACGAGATCTAGTACGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT   CGTACTAG
P7.S703  CAAGCAGAAGACGGCATACGAGATTTCTGCCTGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT   AGGCAGAA
P7.S704  CAAGCAGAAGACGGCATACGAGATGCTCAGGAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT   TCCTGAGC
P7.S705  CAAGCAGAAGACGGCATACGAGATAGGAGTCCGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT   GGACTCCT
P7.S706  CAAGCAGAAGACGGCATACGAGATCATGCCTAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT   TAGGCATG
```

Barcode-swapping blockers (Table S8):

```
Short SBS            CGTGTGCTCTTCCGATCT
Full SBS             CAAGCAGAAGACGGCATACGAGATtcgccttaGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
Decoy DNA strand A   GGTAGAAG/ideoxyU//ideoxyU/AGTAGAATGAAG/ideoxyU//ideoxyU/AGAAGA/ideoxyU//ideoxyU/GTAA/3InvdT/
Decoy DNA strand B   TTACAATC/ideoxyU//ideoxyU/CTAACTTCA/ideoxyU//ideoxyU/CTACTAAC/ideoxyU//ideoxyU/CTACC/3InvdT/
Blocking oligo       CTGTCTCTTATACACATCTCATCATAGAGATCGGAAGAGCACACG/3InvdT/
```

The 96 Tn5 barcodes (Table S5, well 1 = A1 … well 96 = H12, row-major) run from
`GAACCGCG` (A1) to `ATTGTGAA` (H12) 🟢; the full list is in the source file and is not
copied here.

Not given by the paper (from the 10x kit, 🟡 via upstream): the gel-bead oligo
`P5 · <16-nt GEM barcode> · TCGTCGGCAGCGTC` (s5) and the kit's P5-side Sample Index PCR
primer (upstream writes it as bare `illumina.P5`; the 10x kit sequence itself 🔴 here).

### How they interlock — 🟡 (computed with `lib/`)

- **Tn5ME-A = `nextera.ADAPTOR_S5`** (s5 + ME, 33 nt) — the standard Nextera A adaptor.
- **Tn5MErev = `nextera.ME_RC`** (19 nt), 5'-phosphorylated; anneals to the ME of both A and B.
- **Tn5ME-B** (45 nt) = **`illumina.TRUSEQ_READ2[-18:]`** (`CGTGTGCTCTTCCGATCT`) + 8-nt
  barcode + **`nextera.ME`**. So the B end replaces Nextera's s7 with the 3' 18 nt of the
  TruSeq Read 2 primer — the "partial sequence of i7 TruSeq primer" of the methods.
- **Short SBS = `illumina.TRUSEQ_READ2[-18:]`** = the 5' 18 nt of Tn5ME-B exactly. It primes
  on the complement of the B end, while the 10x bead oligo primes on the complement of s5:
  two primers, so exponential PCR in the droplet. Tm 54.6 °C by `chemdraw.tm` (default
  0.5 µM primer, 50 mM Na⁺; the kit buffer is not that).
- **P7.S70x** (66 nt each) = **`illumina.P7`** + 8-nt index + **`illumina.TRUSEQ_READ2`**
  (34 nt). The index is written as the **reverse complement** of the "Barcode" column for
  all six (S701 carries `TCGCCTTA`, Barcode `TAAGGCGA`). In Sample Index PCR the primer
  anneals with only its 3' 18 nt (the SBS part); the 5' 16 nt `GTGACTGGAGTTCAGA`, the
  index and P7 ride along as a tail.
- The six Barcode values `TAAGGCGA … TAGGCATG` are the Illumina Nextera N701–N706 i7
  sequences 🟡 (recalled, not checked: `lib/` has no N7xx set).
- **Full SBS ≡ P7.S701** — identical bases (only the case differs). It was used *inside*
  the droplets in the blocking test, so a library made with it and then sample-indexed with
  S704 carries two i7 barcodes on different molecules — the trick behind Fig. S18 /
  the "~1/3 instead of 1/16" PCR efficiency estimate. 🟢 for the design, 🟡 for the identity.
- **Blocking oligo** (45 nt) = `nextera.ME_RC` + `CATCATAG` + `revcomp(short SBS)`; i.e.
  the reverse complement of a Tn5ME-B whose barcode would be `CTATGATG`, a sequence **not**
  among the 96 Table S5 barcodes. It therefore pairs with any free Tn5ME-B over 37 of 45 nt
  (barcode positions mismatched) and cannot be extended (3' inverted dT).
- **Decoy DNA**: strands A and B (36 nt each) are exact reverse complements when dU is read
  as T; no ME and no primer site, so whatever Tn5 inserts into them cannot be amplified.
  Purpose of the dU 🔴 (not stated; perhaps to allow USER digestion).
- **Tn5 barcodes**: 96 distinct 8-mers, minimum pairwise Hamming distance 4; upstream's
  `Tn5ME-B_barcode.xlsx` is identical to Table S5 in the same well order.

## 3. Step by step (standard txci-ATAC-seq, lung/liver version)

🟢 from the methods ("txci-ATAC-seq using human lung, mouse lung, and mouse liver tissue
samples"; the brain version differs as noted).

1. **Indexed Tn5**: each well of the iTSM plate holds 5 µL of 500 nM Tn5 loaded with
   Tn5ME-A/Tn5MErev and Tn5ME-B(barcode)/Tn5MErev duplexes — every transposome is a
   mixture of A and B, so tagmentation gives A/A, A/B and B/B fragments 🟡
   ([Tn5](../ref/concepts/tn5-tagmentation.md)).
2. **Nuclei**: frozen nuclei thawed, washed (RSB + 0.1 % Tween-20 + 0.1 % BSA), filtered
   (40 µm Flowmi), counted with DAPI; 20,000 nuclei in 7 µL PBSB per well.
3. **Tagment**: + 13 µL TBS (12.5 µL Illumina Tagment DNA Buffer, 0.25 µL 1 % digitonin,
   0.25 µL 10 % Tween-20) into the 5 µL Tn5 → 25 µL; 37 °C 1 h. No stop buffer: chill on
   ice 5 min (an EDTA/spermidine stop was tried in optimisation and found unnecessary).
   Rows G/H mixed 10,000 human + 10,000 mouse nuclei as an in-run barnyard.
4. **Pool and wash**: pool the plate into TMG (10 mM Tris-acetate pH 7.8, 5 mM Mg-acetate,
   10 % glycerol), spin, filter, spin.
5. **Load**: resuspend in loading buffer **+ 5 µM short SBS** (1 µM final in the GEM
   reaction, preprint wording); 100,000 or 200,000 nuclei in 15 µL into the Chromium.
   The loading buffer changed between versions: preprint "1× TB1 + 1× standard storage
   buffer (Illumina)"; published version gives a defined recipe (10 % glycerol, 20 mM NaCl,
   10 mM Tris-HCl pH 7.5, 0.02 mM EDTA, 0.2 mM DTT, 10 mM Tris-acetate pH 7.6, 5 mM
   Mg-acetate, 10 % DMF) 🟢 — i.e. Tn5 storage buffer + TD-type buffer spelled out 🟡.
6. **GEM generation and barcoding** per CG000209 Rev D steps 2–3, but **8 GEM incubation
   cycles** (step 2.5). In each droplet: 72 °C gap fill-in (kit), then bead oligo
   (P5-GEM-s5) and short SBS amplify A/B fragments exponentially 🟡 (cycle profile from the
   kit, 🔴 here).
7. Optional **10 % / 90 % split** of the GEMs (100 µL → 10 + 90 µL) before cleanup; both
   processed in parallel.
8. **Sample Index PCR** (step 4.1): the kit's Single Index N Set A replaced by 2.5 µL of
   **25 µM P7.S70x** (one per 10x lane); **5 cycles** total.
9. **Double-sided size selection** (kit step 4.2); Qubit; 6 % PAGE.

Variants and optimisation runs 🟢:

- **Brain** (Adey lab): TBS + nuclei 20 µL into the iTSM plate, 37 °C 60 min on a
  ThermoMixer (300 rcf); TMG wash; 14 µL nuclei + **1 µL 75 µM short SBS**; **6** GEM
  cycles; 10 %/90 % split with 1 µL / 2 µL of 10 µM i7 TruSeq primer and **8 / 7** SI-PCR
  cycles.
- **Pilot (unmodified 10x)**: 5,000 nuclei/well in 1.25× Tagment DNA buffer (or 10x ATAC
  Buffer B), 55 °C 30 min, EDTA/spermidine stop; 75,000 nuclei/lane; no SBS; 8 SI-PCR
  cycles. Collision rates 40–47 %.
- **Blocking test** (100,000 nuclei/lane, 20,000 nuclei/well, 37 °C 30 min): *Decoy DNA* —
  2.5 µL 50 µM duplex added after tagmentation, 55 °C 10 min; *Blocking oligo* — 2.5 µL
  100 µM into the GEM master mix; *SBS* — 2.5 µL 25 µM **Full SBS** into the master mix.
  SI-PCR stopped at plateau: 4 cycles (SBS) vs 15 (others). SBS won.

## 4. Phased-txci-ATAC-seq (preprint: Fast-txci-ATAC-seq)

🟢 Same chemistry; only the front end changes. Nuclei are diluted in nuclei-freezing buffer
(NFB: 50 mM Tris pH 8.0, 5 mM Mg-acetate, 25 % glycerol, 0.1 mM EDTA, 5 mM DTT, protease
inhibitor) to 3,175/µL; 6.3 µL (20,000 nuclei) per tube, 8 tubes per sample, flash-frozen.
On library day: thaw on ice, + 13.7 µL transposition buffer (12.5 µL 2× Tagment DNA Buffer,
0.7 µL 10× PBS, 0.25 µL 1 % digitonin, 0.25 µL 10 % Tween-20) + 5 µL 500 nM indexed Tn5;
37 °C 60 min. No washing or counting between thaw and tagmentation. 6 samples × 8 wells =
48 Tn5 barcodes; 50,000 or 100,000 nuclei/lane (i7 S705, S706). The rename "Fast" → "Phased"
is the only difference between preprint and published method text found for this section 🟡.

## 5. Final library — 🟡 (computed; matches upstream exactly)

Only A/B fragments carry both P5 (via the bead) and P7 (via the SBS handle). A/A gets P5 at
both ends, B/B the SBS/i7 at both ends; neither clusters, and upstream adds that both are
suppressed in the droplet by their inverted terminal repeats (semi-suppressive PCR) 🟡.

Top strand, 5'→3' (171 nt + insert):

| # | Segment | Length | Sequence / origin |
|---|---|---|---|
| 1 | P5 | 29 | `illumina.P5` (bead oligo) |
| 2 | GEM barcode | 16 | bead |
| 3 | s5 | 14 | `nextera.S5` (bead oligo 3' end = Tn5ME-A 5' end) |
| 4 | ME | 19 | `nextera.ME` (Tn5ME-A) |
| 5 | insert | variable | genomic DNA |
| 6 | ME' | 19 | `nextera.ME_RC` (from the gap-filled Tn5MErev side) |
| 7 | Tn5 barcode' | 8 | reverse complement of the Table S5 barcode |
| 8 | TruSeq Read 2' | 34 | `revcomp(illumina.TRUSEQ_READ2)`; 18 nt from Tn5ME-B, 16 nt from the P7.S70x tail |
| 9 | i7 | 8 | the Table S7 "Barcode" column as written |
| 10 | P7' | 24 | `illumina.P7_RC` |

Assembled from `lib/` this is character-for-character upstream's "final library
structure" (with N for barcodes).

## 6. Sequencing

🟢 Final runs (lung/liver): NextSeq 550 High Output, **Read 1 51, i7 10, i5 16, Read 2 78**
(only 8 i7 cycles needed; 10 run to match a co-pooled bulk ATAC library). Pilot: 50/8/16/77.
Blocking test: 50/10/16/92. Brain: NextSeq 500 Mid or NovaSeq 6000 S4 (cycles not given).

| Read | Primer 🟡 | Reads |
|---|---|---|
| Read 1 (51) | Nextera Read 1 = `nextera.READ1_PRIMER` | genomic DNA from the A end |
| i7 (8) | TruSeq index 1 = `illumina.INDEX1_PRIMER` | i7 = Table S7 Barcode column (computed) |
| i5 (16) | `nextera.INDEX2_PRIMER` (upstream's "GEM barcode sequencing primer") | 10x GEM barcode |
| Read 2 (78) | TruSeq Read 2 = `illumina.TRUSEQ_READ2` | 8 nt Tn5 barcode (as written in Table S5, computed) + 19 nt ME + 51 nt gDNA |

🟢 The analysis takes the first 8 bp of Read 2 as the Tn5 barcode and trims 27 bp (8 + 19 ME),
correcting barcodes within edit distance 2. Bases 9–27 of Read 2 are the fixed ME, so the
authors **spike 5 % bulk ATAC** for base diversity and recommend PhiX or **dark cycles 9–27**
when sequenced alone 🟢. Read primers are the instrument's standard mixes — "to avoid
sequencing with a custom recipe" is the stated reason for the B-adaptor design 🟢; which
standard primer mix covers both Nextera Read 1 and TruSeq Read 2 is 🟡.

## 7. Upstream (scg_lib_structs) checked against the paper

| Upstream claim | Paper | Verdict |
|---|---|---|
| Tn5ME-A, Tn5ME-B, Tn5MErev sequences | Table S4 | agree; upstream writes `/phos/`, paper `[phos]` |
| Tn5 barcode plate | Table S5 | identical, same well order (computed) |
| Custom i7 TruSeq primer `P7 · [8-bp i7] · TruSeq R2` | Table S7 | agree; upstream does not say the i7 is written reverse-complemented in the oligo |
| "SBS Short Primer" added for exponential in-droplet PCR | short SBS in final protocol; **Full SBS** in the blocking test | agree for the final protocol; upstream omits the Full SBS, blocking oligo and decoy DNA |
| Gap fill 72 °C 5 min | "per 10x CG000209" only | 🟡 (kit) |
| Final structure and read lengths 51/8/16/78 | 51/10/16/78 (8 i7 cycles suffice) | agree |
| R2 = 8 bp Tn5 barcode + 19 ME + 51 gDNA | first 8 bp barcode, 27 trimmed | agree |
| "custom sequencing primers are not needed" | stated rationale | agree |
| Bead oligo, sequencing primers | not in paper | 🟡 from upstream / 10x |
| Tn5 stop, GEM cycle numbers, Phased variant, 10 %/90 % split | in paper | not on the upstream page |

## 8. Open questions

- 🔴 The protocols.io protocol was not read; volumes there may differ from the paper.
- 🔴 Why the Decoy DNA carries pairs of deoxyuridines.
- 🔴 The blocking oligo's barcode `CTATGATG` — a deliberate non-barcode, or a leftover
  from another plate?
- 🔴 The 10x kit's in-droplet cycling profile and the exact P5-side SI-PCR primer.
- 🟡 Why exponential PCR helps: the paper argues it outcompetes swapping; the mechanism
  of swapping (free Tn5ME-B priming on other fragments in the droplet vs. something else)
  is only sketched in Fig. S1c.
- 🟡 SBS shifts the fragment-size distribution toward short fragments (main text; fewer
  GEM cycles partly restore it — Fig. S2g, 3 vs 8 cycles).

## 9. How this note was made (tool evaluation)

`tools/get_sources.py` fetched only the preprint PDF and the upstream page; the Genome
Biology version and its Additional files were found and fetched by hand (Europe PMC
XML and Springer static-content URLs), as was upstream's barcode xlsx. The preprint PDF
happens to print Tables S4, S7, S8 at the end, so the oligos were available either way;
Table S5 (barcodes) is only in the separate files. `tools/scrape_primers.py` found all
linker, i7 and blocker oligos with `illumina.P7` / `TRUSEQ_READ2` / `nextera.ME` already
recognised. `--find` does not match a query containing `N` (Tn5ME-B had to be checked in
two halves) and missed the lowercase-barcode Full SBS rows of Table S8 (it reported only the
identical P7.S701 rows).
