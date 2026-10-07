# sci-RNA-seq and sci-RNA-seq3: split-pool combinatorial indexing of fixed cells or nuclei

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published (or not obtainable here). Relationships marked 🟡 *(computed)* were worked out
> with `lib/` while writing this note.
> Claims taken only from the upstream scg_lib_structs page are 🟡 *(upstream)*: a
> secondary source.

**sci-RNA-seq**: Cao J, Packer JS, Ramani V, Cusanovich DA, Huynh C, Daza R, Qiu X,
Lee C, Furlan SN, Steemers FJ, Adey A, Waterston RH, Trapnell C, Shendure J.
"Comprehensive single-cell transcriptional profiling of a multicellular organism."
*Science* 357:661–667 (2017). doi:[10.1126/science.aam8940](https://doi.org/10.1126/science.aam8940),
PMC5894354. Preprint: *bioRxiv* doi:[10.1101/104844](https://doi.org/10.1101/104844)
(v1, 2 Feb 2017; CC-BY-NC-ND 4.0).

**sci-RNA-seq3**: Cao J, Spielmann M, Qiu X, Huang X, Ibrahim DM, Hill AJ, Zhang F,
Mundlos S, Christiansen L, Steemers FJ, Trapnell C, Shendure J. "The single-cell
transcriptional landscape of mammalian organogenesis." *Nature* 566:496–502 (2019).
doi:[10.1038/s41586-019-0969-x](https://doi.org/10.1038/s41586-019-0969-x), PMC6434952.

Sources read (fetched by `tools/get_sources.py`, into
`$CHEM_DATA/sources/sci-rna-seq-family__10.1126+science.aam8940/`, never committed):

| File | What | Used for |
|---|---|---|
| `104844_v1.full.pdf` (.txt) | sci-RNA-seq preprint, with Materials & Methods | **all sci-RNA-seq chemistry**: RT primer, PCR primers, conditions, read lengths |
| `science.aam8940_PMC5894354.html` (.txt) | Science author manuscript on PMC | workflow summary, Fig. 1 legend, the three-level (indexed Tn5) variant; it has **no methods** (they are in the gated supplement) |
| `s41586-019-0969-x_PMC6434952.xml` (.txt) | sci-RNA-seq3, Europe PMC JATS full text | **all sci-RNA-seq3 methods** |
| `s41586-019-0969-x_41586_2019_969_MOESM3_ESM.xlsx` (.txt) | Nature Supplementary Tables 1–11 | **Table S11: every sci-RNA-seq3 oligo** (384 RT, 384 ligation, 96 P5, 96 P7) |
| `s41586-019-0969-x_41586_2019_969_MOESM1_ESM.pdf`, `..._MOESM2_ESM.pdf` | Supplementary Note 1 (Monocle), Reporting Summary | not chemistry |
| `upstream_sci-RNA-seq_family.html` (.txt) | Teichlab scg_lib_structs page | second source, checked against the above (§7) |

Could not be fetched (get by hand and drop in the same directory):

| What | URL | Why it matters |
|---|---|---|
| Science supplement (NIHMS956226, Supplementary Materials, 2.2 MB PDF) | https://pmc.ncbi.nlm.nih.gov/articles/instance/5894354/bin/NIHMS956226-supplement-supplement_1.pdf | the **published** sci-RNA-seq methods, the RT barcode list, the i5/i7 list, the **indexed Tn5 adaptors** of the three-level variant, Table S1 (experiments) 🔴 |
| Science version of record, with supplement | https://www.science.org/doi/10.1126/science.aam8940 | same; Cloudflare / paywall |

`s41586-019-0969-x_PMC6434952.html` is a reCAPTCHA page, not the article, and
`s41586-019-0969-x_PMC6434952_supplementary.zip` is truncated (no central directory); the
JATS XML and the Springer supplements replace both.

---

## 1. What it is

No compartments, no beads: **the cell (or nucleus) itself is the reaction vessel**, and its
identity is the combination of well-specific barcodes it collected while being split,
pooled and split again. 🟡 (paraphrase: the Science main text contrasts it with methods
that isolate cells in physical compartments; workflow per Fig. 1A)

- **sci-RNA-seq** (2017), two levels: methanol-fixed cells or nuclei are split into
  96/384 wells for **in situ RT with a barcoded oligo-dT** (barcode 1 + UMI), pooled,
  **FACS-sorted** 10–100 per well into a second plate, and given second-strand synthesis,
  Tn5 tagmentation and **indexed PCR** (barcode 2 = i5 + i7 pair). 🟢 An optional third
  level puts a barcode on **indexed Tn5** between RT and PCR (16 × 6 × 16 demonstrated). 🟢
  (Science main text; adaptor sequences 🔴, in the gated supplement)
- **sci-RNA-seq3** (2019), three levels: nuclei straight from tissue, PFA-fixed; RT
  barcode, then a **barcoded hairpin adaptor ligated to the 5' end of the RT primer**
  (replacing indexed Tn5 as the third level), then second-strand synthesis, tagmentation
  with **Tn5 loaded only with the N7 (s7) adaptor**, **USER** cleavage, indexed PCR. Sorting
  is replaced by dilution. 384 × 384 × 768 ≈ 1.1 × 10⁸ combinations. 🟢 (Nature main text
  and Methods)

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | **Combinatorial indexing of whole fixed cells**: the first barcode is put on in situ, the cell stays intact through pooling | sci-ATAC-seq, sci-CAR, SPLiT-seq / SHARE-seq (ligation rounds), sci-Plex |
| 2 | **Read 1 handle on the RT primer, s7 from Tn5 at the other end**: only the 3' end of each cDNA, between the RT primer and the nearest s7 insertion, amplifies | the same 3'-tag logic as CEL-seq / Drop-seq, built from TruSeq + Nextera halves |
| 3 | (sci3) A **6-nt splint-ligated hairpin** adds a barcode to the 5' end of the RT primer *inside* the nucleus, with a **deoxyuridine** to cut the hairpin open again before PCR | the dU-hairpin trick of the NEBNext adaptor (`illumina.NEBNEXT_HAIRPIN`) |
| 4 | (sci3) **N7-only Tn5**: every insertion carries s7, so none is wasted on s5 ends | "Tn5 transposomes loaded only with N7 adaptor increased UMI counts by over 50%" (Extended Data Fig. 1b) 🟢 |

Concepts: [reverse transcription](../ref/concepts/reverse-transcription.md) (oligo-dT
priming, RNase H / Pol I second strand),
[Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md). No template switching: the 5' end
of the transcript is never captured.

## 2. Oligos

### sci-RNA-seq (2017) 🟢

Verbatim from the preprint Materials & Methods (`104844_v1.full.pdf.txt` lines 899–900
and 922–931; the PDF wraps them, re-joined here). Brackets are the paper's placeholders.

```
anchored oligo-dT (RT)  5′-ACGACGCTCTTCCGATCTNNNNNNNN[10bp index]TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTVN-3′   (IDT)
P5 primer               5′-AATGATACGGCGACCACCGAGATCTACAC[i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT-3′
P7 primer               5′-CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGG-3′
Tn5                     Nextera DNA Sample Preparation kit, TDE1 (standard s5/s7 loading)
```

"N" is any base, "V" is A, C or G. 🟢 The individual RT barcodes, the i5 / i7 sequences
and the indexed-Tn5 adaptors of the three-level variant are **not** in the preprint or
the PMC main text 🔴 (gated Science supplement, see table above). The same RT primer is
used for the bulk "Tn5-RNA-seq" comparison library (preprint line 1127). 🟢

### sci-RNA-seq3 (2019) 🟢

Verbatim from Nature Supplementary Table 11 (`..._MOESM3_ESM.xlsx.txt`, sheet "Table S11",
"RT, ligation and PCR primer sequences for sci-RNA-seq3"). IDT notation as given:
`/5Phos/` = 5' phosphate, `/ideoxyU/` = internal 2'-deoxyuridine. The P7 table writes
the i7 in **lower case**.

```
sc_ligation_RT_1     /5Phos/CAGAGCNNNNNNNNTCCTACCAGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
sc_ligation_RT_384   /5Phos/CAGAGCNNNNNNNNTTGGAGTCGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
   ... 384 RT primers (sc_ligation_RT_1..384)

sc_ligation_1        GCTCTGACAATCAAGT/ideoxyU/ACGACGCTCTTCCGATCTACTTGATTGT
sc_ligation_382      GCTCTGATGAGCCAA/ideoxyU/ACGACGCTCTTCCGATCTTTGGCTCAT
   ... 384 hairpin ligation adaptors (sc_ligation_1..384)

P5-1                 AATGATACGGCGACCACCGAGATCTACACCTCCATCGAGACACTCTTTCCCTACACGACGCTCTTCCGATCT
   ... P5-1..96
P7-1                 CAAGCAGAAGACGGCATACGAGATccgaatccgaGTCTCGTGGGCTCGG
   ... P7-1..96

Tn5                  "i7 only TDE1 enzyme (62.5 nM, Illumina)"
```

The Methods text gives the same oligos as templates 🟢 (`..._PMC6434952.xml.txt`, line 47–51):
RT `5′-/5Phos/CAGAGCNNNNNNNN[10bp barcode]TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT-3′`; ligation
adaptor `5’-GCTCTG[9 bp or 10 bp barcode A]/dideoxyU/ACGACGCTCTTCCGATCT[reverse complement
of barcode A]-3’`; P5 `5′-AATGATACGGCGACCACCGAGATCTACAC[i5]ACACTCTTTCCCTACACGACGCTCTTCCGATCT-3′`;
P7 `5′-CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGG-3′`.

Note the modification name: the text says **"dideoxyU"**, the oligo table says
**`/ideoxyU/`** (IDT's internal deoxyuridine). Only a deoxyuridine can sit internally in
an oligo and be excised by USER; a dideoxy nucleotide has no 3'-OH and could only be
terminal. The table is right, the text is a slip. 🟡

### How they interlock — 🟡 (computed with `lib/`, over all 960 rows of Table S11)

**One 18-nt handle, `ACGACGCTCTTCCGATCT`, is the bridge in both versions.** It is
`illumina.TRUSEQ_READ1[-18:]`, the 3' end of the 33-nt TruSeq Read 1 primer site.

- sci-RNA-seq: the RT primer starts with it. The P5 primer =
  **`illumina.P5` + [i5] + `illumina.TRUSEQ_READ1`** (33 nt), so the primer's last 18 nt
  are the handle itself and it primes on the handle's complement in the second strand. The
  first PCR cycle adds the other 15 nt of Read 1 and P5.
- sci-RNA-seq3: the handle has moved off the RT primer into the **loop of the hairpin
  adaptor**; the P5 primers are unchanged. All 96: **`illumina.P5` + 10-nt i5 +
  `illumina.TRUSEQ_READ1`**, 72 nt.
- P7 primer (both versions) = **`illumina.P7` + [i7] + `nextera.S7`**: a Nextera N7xx
  primer with a 10-nt index. Its 3' end is s7 only (15 nt, no ME), so it primes on
  the s7 complement left by tagmentation and gap fill. All 96 sci3 P7 primers fit this,
  49 nt, i7 10 nt. Primer-arm Tm (`chemdraw.tm`, 50 mM Na⁺): s7 54.4 °C, the 18-nt
  handle 55.9 °C, full TruSeq Read 1 66.9 °C; annealing in both papers is 66 °C, which
  the P7 arm only reaches once the whole primer is incorporated after cycle 1.

sci-RNA-seq3 RT primers (all 384): **`/5Phos/` + `CAGAGC` + N₈ (UMI) + 10-nt RT barcode +
T₃₀**, 54 nt, **no VN anchor** (the 2017 primer has one), 384 distinct barcodes, minimum
pairwise Hamming distance 3. 87 of the barcodes end in T, so the oligo-dT run is
effectively 31–32 T for those; the table and the text agree on T₃₀ after the barcode.

sci-RNA-seq3 hairpin adaptors (all 384): **`GCTCTG` + bcA + `/ideoxyU/` + handle +
revcomp(bcA)**. In all 384 the 3' arm is exactly the reverse complement of bcA. bcA is
9 nt in 238 adaptors (43-mer) and 10 nt in 146 (45-mer), and **every bcA starts with A**,
so every adaptor ends in T. The 384 read-side barcodes (revcomp(bcA)) are distinct.

How the hairpin finds the RT primer:

- `GCTCTG` is the reverse complement of **`CAGAGC`**, the 5' end of every RT primer.
- Folded, the adaptor is a stem of bcA : revcomp(bcA) closed by a 19-nt loop (dU +
  handle), with **`GCTCTG` as a 6-nt 5' overhang**.
- The overhang pairs with the 5'-phosphorylated `CAGAGC` of the first-strand cDNA. That
  puts the RT primer's 5'-phosphate right against the adaptor's 3'-OH (revcomp(bcA)), a
  nick that T4 DNA ligase seals. The product is a 15–16 bp stem
  (GCTCTG-bcA : revcomp(bcA)-CAGAGC) closing over the dU + handle loop, continuous with
  the cDNA.

That is why the RT primer is phosphorylated and the 2017 handle is gone from it: in sci3
the RT primer's 5' end is a ligation splint site, not a PCR handle.

## 3. Step by step

### sci-RNA-seq (preprint methods) 🟢 unless marked

1. **Fix / isolate**: cells in ice-cold 100 % methanol, −20 °C 10 min, DEPC and
   SUPERase-In / BSA washes; or nuclei by IGEPAL lysis. 5000 cells or nuclei per µL.
2. **In situ RT, barcode 1**: per well 1,000–10,000 cells (2 µL) + 1 µL 25 µM barcoded
   oligo-dT + dNTPs, 55 °C 5 min, ice; SuperScript IV mix; **55 °C 10 min**; stopped with
   2× stop solution (40 mM EDTA, 1 mM spermidine). First strand (🟡):
   `5'- handle · UMI · RT barcode · T30 · V · N · <cDNA, antisense> -3'`.
3. **Pool, stain DAPI, FACS** 10–100 cells per well (an example range; per-experiment numbers are in Table S1 🔴)
   into 5 µL EB; DAPI gating removes doublets.
4. **Second strand**: NEB mRNA Second Strand Synthesis module (RNase H / *E. coli* Pol I /
   ligase), 16 °C 150 min, 75 °C 20 min. 🟡 (module contents from NEB, not the paper)
5. **Tagment**: + **5 ng human genomic DNA as carrier** (against over-tagmentation and
   loss), Nextera TD buffer, 0.5 µL TDE1, 55 °C 5 min; Zymo DNA binding buffer, AMPure
   (36 µL), elute 16 µL. The methods state that carrier DNA is not appreciably amplified
   because the PCR primers are specific to the RT products 🟢; the mechanism (carrier has
   s5/s7 ends only, and the P5 primer needs the handle) is 🟡.
6. **Indexed PCR, barcode 2**: P5[i5] + P7[i7], NEBNext HiFi 2× master mix;
   **75 °C 3 min** (gap fill of the 9-nt Tn5 gap 🟡), 98 °C 30 s, **18–22 ×** (98 °C 10 s,
   66 °C 30 s, 72 °C 1 min), 72 °C 5 min. Pool, AMPure 0.8×, Qubit, 6 % TBE-PAGE.

   Which fragments amplify (🟡 *(computed)*, same as the upstream drawing): only the one
   from the **RT-primer end to the nearest s7 insertion**. Handle-to-s5 fragments, s5–s5,
   s7–s7 and s5–s7 internal fragments (and all carrier DNA) lack one of the two primer
   sites, since the P5 primer needs the TruSeq handle and does not bind s5.
7. **Three-level variant** (Science only): after RT and second-strand synthesis, pool,
   split to 6 wells for **indexed Tn5**, pool, sort to 16 wells for PCR. 🟢 Adaptor
   sequences, where the Tn5 barcode lands in the read, and the sequencing recipe are in
   the Science supplement 🔴. The paper notes fewer UMIs per cell than two-level at equal
   depth, possibly from lower efficiency of indexed Tn5. 🟢

### sci-RNA-seq3 (Nature methods) 🟢 unless marked

The published Methods swap two phrases: the paragraph that **sets up the ligation** begins
"After ligation reaction…", and the one that **follows the ligation** begins "After RT
reaction…". The order of the plates (4 RT plates, then 4 ligation plates, then 8 plates for
second strand / tagmentation / PCR) and the 384 × 384 × 768 arithmetic fix the order as
RT → ligation → second strand. 🟡 (Upstream has the same order.)

1. **Nuclei**: extracted from minced embryo in cell lysis buffer, fixed in 4 %
   paraformaldehyde 15 min on ice, flash frozen. Thawed: 0.2 % Triton X-100 3 min,
   **brief sonication** (12 s, low), Flowmi filtration.
2. **In situ RT, barcode 1** (4 × 96 wells): 80,000 nuclei (16 µL) + 8 µL 25 µM RT primer
   + dNTPs, 55 °C 5 min, ice; SuperScript IV mix (14 µL); **gradient 4, 10, 20, 30, 40,
   50 °C × 2 min each, then 55 °C 10 min**. First strand (🟡):
   `5'-P · CAGAGC · UMI · RT barcode · T30 · <cDNA, antisense> -3'`.
3. **Pool**, spin 500 × g, redistribute to 4 × 96 wells.
4. **Hairpin ligation, barcode 2**: per well 4 µL T4 ligase buffer, 2 µL T4 DNA ligase,
   4 µL 5 M betaine, 6 µL nuclei, **8 µL 100 µM barcoded hairpin adaptor**, 16 µL 40 %
   PEG 8000; **16 °C 3 h**. Product (🟡 *(computed)*, written as the continuous ligated
   strand):
   `5'- GCTCTG · bcA · dU · handle(18) · revcomp(bcA) · CAGAGC · UMI · RT bc · T30 · <cDNA> -3'`,
   folded back on itself between GCTCTG-bcA and revcomp(bcA)-CAGAGC.
5. **Pool, dilute, filter**: 60 µL dilution buffer, pool, 600 × g, Flowmi twice, count,
   **2,500 nuclei per well** into 8 × 96 wells (5 µL wash buffer + 5 µL EB). No FACS.
6. **Second strand**: 1.33 µL NEB mRNA Second Strand buffer + 0.66 µL enzyme, 16 °C 180 min.
   The second strand runs into the hairpin and, by strand displacement, copies it
   to the end, so the second strand carries the complement of the whole hairpin with an A
   opposite the dU. 🟡 (inferred; upstream draws the same; the paper does not describe it)
7. **Tagment**: 11 µL Nextera TD buffer + 1 µL **"i7 only" TDE1** (62.5 nM), 55 °C 5 min;
   24 µL Zymo DNA binding buffer, 5 min; AMPure 1.5×. Every insertion carries s7-ME.
8. **USER, on the beads**: elute in 8 µL water + 1 µL 10× USER buffer + 1 µL USER enzyme,
   37 °C 15 min; +6.5 µL EB; remove beads. USER excises the dU and nicks the first strand
   there, so the 5' `GCTCTG · bcA` piece (15–16 nt) is cut off the strand that carries the
   handle. 🟡 Why (most likely so that the handle is no longer locked in a hairpin loop that
   can snap back, keeping the P5 primer site accessible) is not stated 🔴.
9. **Indexed PCR, barcode 3**: 16 µL + 2 µL 10 µM P5-n + 2 µL 10 µM P7-n + 20 µL NEBNext
   HiFi 2×; **72 °C 5 min** (gap fill 🟡), 98 °C 30 s, **12–14 ×** (98 °C 10 s, 66 °C 30 s,
   72 °C 1 min), 72 °C 5 min. Pool, AMPure 0.8×, Qubit, 6 % TBE-PAGE.

   "768 barcodes introduced by PCR" 🟢 = 8 plates × 96 wells, each well a distinct P5/P7
   pair; the 96 + 96 primers of Table S11 allow up to 9,216 pairs. Which pairs were used is
   not given. 🔴

## 4. Final libraries — 🟡 (assembled from the oligos above; agrees with upstream)

Top strand, 5'→3'. The cDNA is the first (antisense) strand, read from the poly(T) inward.

**sci-RNA-seq** (two-level):

```
5'- P5 · i5(10) · TRUSEQ_READ1 [last 18 nt = RT handle] · UMI(8) · RT barcode(10) · T30 · V · <cDNA 3' end, antisense> · ME_RC · S7_RC · i7'(10) · P7_RC -3'
```

**sci-RNA-seq3** (three-level):

```
5'- P5 · i5(10) · TRUSEQ_READ1 [last 18 nt = hairpin loop handle] · revcomp(bcA)(9|10) · CAGAGC · UMI(8) · RT barcode(10) · T30 · <cDNA 3' end, antisense> · ME_RC · S7_RC · i7'(10) · P7_RC -3'
```

`ME_RC` = `nextera.ME_RC`, `S7_RC` = `nextera.S7_RC`, `P7_RC` = `illumina.P7_RC`; i7' is
the reverse complement of the i7 as written in the P7 primer. The ME and the part of s7
beyond the primer come from the Tn5 adaptor; the 10-nt i7 and P7 from the PCR primer.

## 5. Read layout / sequencing

| | sci-RNA-seq 🟢 | sci-RNA-seq3 🟢 |
|---|---|---|
| Instrument | NextSeq 500, V2 75-cycle kit | NovaSeq |
| Read 1 | 18 cycles | 34 cycles |
| Index 1 (i7) | 10 | 10 |
| Index 2 (i5) | 10 | 10 |
| Read 2 | 52 | 52 |

What each read covers, 🟡 *(computed from §4)*:

- **Read 1** is primed by the standard TruSeq Read 1 primer (`illumina.TRUSEQ_READ1`,
  which the library carries whole) and reads cell barcodes only:
  - sci-RNA-seq: UMI(8) + RT barcode(10) = 18, exactly the 18 cycles; the preprint Fig. 1b
    legend says the same 🟢.
  - sci-RNA-seq3: revcomp(bcA)(10) + `CAGAGC` + UMI(8) + RT barcode(10) = 34. With a 9-nt
    bcA it is 33 and the 34th cycle reads the first T, so the RT barcode does not sit at a
    fixed position: the parser must find the ligation barcode first (the paper
    demultiplexes on "RT index and ligation index", ED < 2 including indels 🟢).
- **Read 2** is primed by the Nextera Read 2 primer (`nextera.READ2_PRIMER` = s7 + ME,
  present via the Tn5 adaptor) and reads 52 nt of the transcript in the **sense**
  orientation, starting at the Tn5 insertion site, which is the UMI-collapsing coordinate
  ("read 2 end-coordinate" 🟢).
- **Index 1** reads i7 through `nextera.INDEX1_PRIMER` (ME_RC + S7_RC), giving the
  **reverse complement** of the i7 written in the P7 primer (P7-1 `ccgaatccga` reads
  `TCGGATTCGG`).
- **Index 2** reads i5 next to TruSeq Read 1 (upstream gives the reverse-complement
  workflow primer `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` = `illumina.INDEX2_PRIMER_RC`);
  orientation depends on the instrument workflow (NextSeq 500 reverse complement;
  NovaSeq 6000 v1.0 forward).
- Cell = RT barcode + PCR i5/i7 pair (sci-RNA-seq), plus the ligation barcode (sci3).
  None of the read primers is custom; the papers name none. 🟡

## 6. What upstream (scg_lib_structs) says, checked

| Upstream claim | Primary source | Verdict |
|---|---|---|
| sci-RNA-seq RT primer `ACGACGCTCTTCCGATCT[8-bp UMI][10-bp RT barcode]T30VN` | preprint, line 899 | **agrees** |
| P5 / P7 primers, s7, carrier gDNA, read lengths 18/10/10/52 | preprint methods | **agrees** |
| i5 and i7 drawn as 10 N | 10 index cycles (both); Table S11 (sci3) | **agrees** |
| sci3 RT primer `/Phos/CAGAGC[8-bp UMI][10-bp RT barcode]` T30 **VN** | Methods text and all 384 rows of Table S11: T30, **no VN** | **disagrees**: upstream adds a VN the sci3 oligos do not have |
| hairpin `GCTCTG[reverse complement of barcode A]/ddU/ACGACGCTCTTCCGATCT[9-bp or 10-bp barcode A]` | text: `GCTCTG[barcode A]/dideoxyU/…[reverse complement of barcode A]`; table: `/ideoxyU/` | same molecule with the barcode naming swapped; both upstream (`ddU`) and the paper text ("dideoxyU") misname the internal **deoxy**uridine of the oligo table |
| 384 RT primers, 384 hairpins, 96 P5, 96 P7 | Table S11 | **agrees** |
| sci3 Read 1 34 cycles = hairpin barcode + `GTCTCG` + UMI + RT barcode | the bases read are `CAGAGC` (top strand); `GTCTCG` is the bottom strand written 3'→5' | wording only |
| sci3 Tn5 "homodimer with s7-ME oligo" | "i7 only TDE1" | **agrees** |
| order RT → ligation → second strand → tagment → USER → PCR | Methods (with the swapped phrases) | **agrees** |
| sci-RNA-seq demonstrated as two-level; three-level by indexed Tn5 "can be used" | Science main text | agrees; upstream does not draw the three-level structure either |

## 7. Open questions

- 🔴 The sci-RNA-seq RT barcodes, i5/i7 lists and the **indexed Tn5 adaptors** (three-level
  variant): in the gated Science supplement. Is the Tn5 barcode read in Index 1 or
  inside Read 2?
- 🔴 Why sci3 drops the **VN anchor** of the 2017 RT primer.
- 🔴 Why every hairpin barcode starts with **A** (so the adaptor's 3' end, the ligated
  base, is always T) and why 9- and 10-nt barcodes are mixed. A staggered length would
  shift Read 1 phasing on low-diversity cycles; not stated.
- 🔴 The purpose of **USER** cleavage (step 8): inferred above, not explained in the paper.
- 🔴 Which 768 of the 9,216 possible P5/P7 pairs were used, and whether i5 is read
  forward or reverse complement on the NovaSeq run.
- 🟡 Methods phrase swap ("After ligation" / "After RT") read as an editing slip; the plate
  arithmetic supports RT → ligation.

## 8. How this note was made (tool evaluation)

`tools/get_sources.py` got the bioRxiv preprint (which, unlike the Science PMC
manuscript, carries Materials & Methods), the Europe PMC XML and all three Springer
supplements for sci-RNA-seq3, and listed the Science supplement as `(manual)`. The PMC
HTML for the Nature paper came back as a reCAPTCHA page and the Europe PMC supplementary
zip was truncated; neither mattered. `tools/scrape_primers.py --max-hits 80` put Table S11
first and parsed `/5Phos/` and `/ideoxyU/` correctly, but its 80 hits for that file were
RT, ligation and P5 rows only: **none of the 96 P7 primers** (lower-case i7) was listed,
and in the text it found the polyT but not the `/5Phos/CAGAGCNNNNNNNN` 5' part of the
sci3 RT template. All 960 rows of the sheet were therefore parsed by hand with `lib/`.
`--find` could not confirm sequences containing an `N` run followed by bases, or with a
lower-case index (`P7-…`); those were confirmed with a plain text search of the `.txt`
twin (Table S11 rows `sc_ligation_RT_1`, `sc_ligation_RT_384`, `P7-1`).
