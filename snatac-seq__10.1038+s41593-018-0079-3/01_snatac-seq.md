# snATAC-seq: two-level combinatorial ATAC-seq on nuclei from frozen tissue

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published / not available. Relationships marked 🟡 *(computed)* were worked out with
> `lib/` while writing this note; they are not yet asserted in a self-test, because this
> protocol has no `tools/` module yet (status `notes`).

**snATAC-seq**: Preissl S, Fang R, Huang H, Zhao Y, Raviram R, Gorkin DU, Zhang Y, Sos BC,
Afzal V, Dickel DE, Kuan S, Visel A, Pennacchio LA, Zhang K, Ren B. "Single-nucleus
analysis of accessible chromatin in developing mouse forebrain reveals cell-type-specific
transcriptional regulation." *Nature Neuroscience* 21:432–439 (2018).
doi:[10.1038/s41593-018-0079-3](https://doi.org/10.1038/s41593-018-0079-3), PMID 29434377,
PMC5862073. Analysis pipeline: github.com/r3fang/snATAC (also shipped as Supplementary
Software).

Related papers (from `catalogue/scg_lib_structs.tsv`):

- **sci-ATAC-seq**: Cusanovich DA *et al.* *Science* 348:910 (2015),
  doi:[10.1126/science.aab1601](https://doi.org/10.1126/science.aab1601). The parent
  method; snATAC-seq calls itself "combinatorial ATAC-seq … with modifications". See
  [the sci-ATAC-seq note](../sci-atac-seq__10.1126+science.aab1601/01_sci-atac-seq.md).
- **CPT-seq / Amini 2014** (*Nat. Genet.* 46:1343, doi:10.1038/ng.3119): the source of
  the barcoded transposon and sequencing-primer sequences. The snATAC-seq oligo table
  says so in a footnote 🟢: all oligos except the spike-ins are "from Amini et al., Nat.
  Genet., 2014".

Sources read (fetched by `tools/get_sources.py`, into
`_data/sources/snatac-seq__10.1038+s41593-018-0079-3/`, never committed):

| File | What | Used for |
|---|---|---|
| `s41593-018-0079-3_PMC5862073.html.txt` | full text incl. Online Methods (author manuscript on PMC) | **all steps, conditions, read layout** |
| `s41593-018-0079-3_41593_2018_79_MOESM5_ESM.xlsx.txt` | Supplementary Table 5 (sheets "PCR Primer", "Tagmentation Oligos") | **every oligo** |
| `s41593-018-0079-3_41593_2018_79_MOESM1_ESM.pdf.txt` | Supplementary Figures 1–4 legends | optimisation of the protocol (Supp. Fig. 1) |
| `s41593-018-0079-3_41593_2018_79_MOESM2_ESM.pdf.txt` | Reporting summary | NeuN-sort sample prep, Tn5 source |
| `s41593-018-0079-3_41593_2018_79_MOESM3_ESM.pdf.txt` | Supplementary Tables 1–3 legends | barcode "Set_1 / Set_2" per replicate, comparison table |
| `moesm6_unzipped/snATAC-master/` | Supplementary Software (pipeline) | `barcodes/{i5,i7,r5,r7}_ATAC`: barcode orientation as read |
| `upstream_snATAC-seq.html.txt` | scg_lib_structs page | short text; checked below |
| `science.aab1601_PMC4836442.html.txt` | sci-ATAC-seq 2015 paper | context only (covered in its own note) |
| `…MOESM4_ESM.xlsx.txt` | Supplementary Table 4 (peaks, motifs) | not chemistry, not read |

Not fetched (`(manual)` in the MANIFEST, PMC download gate). They are the NIHMS copies of
the same supplements, so probably nothing is missing:

- https://pmc.ncbi.nlm.nih.gov/articles/instance/5862073/bin/NIHMS950173-supplement-Supplemental_table_5.xlsx
  (the oligo table again, same as MOESM5 presumably)
- https://pmc.ncbi.nlm.nih.gov/articles/instance/5862073/bin/NIHMS950173-supplement-supplemental_text_and_figures.pdf
- https://pmc.ncbi.nlm.nih.gov/articles/instance/4836442/bin/NIHMS776051-supplement-Supplementary_Material.pdf
  (sci-ATAC-seq 2015 supplement: only needed for the parent method)

---

## 1. What it is, and what is new

Two rounds of barcoding on intact nuclei, the sci-ATAC-seq scheme 🟢:

1. **Round 1, in the Tn5.** 96 wells of 4,500 nuclei each; each well gets a Tn5 loaded
   with one of 8 "A" (p5) transposons and one of 12 "B" (p7) transposons, so the well is
   encoded by a pair of 8-nt barcodes inserted at both ends of every fragment.
2. **Round 2, in the PCR.** Nuclei are pooled, and 25 are sorted into each of 384 wells;
   a Nextera-style i5 + i7 primer pair per well adds two more 8-nt indices.

A nucleus is identified by four barcodes: p5 + p7 (tagmentation well) and i5 + i7 (PCR
well). There is no UMI.

What snATAC-seq changes relative to sci-ATAC-seq 🟢 (Methods; Supp. Fig. 1):

| | Change | Why (as given) |
|---|---|---|
| 1 | **Frozen tissue**: lysis in NPB (5 % BSA, 0.2 % IGEPAL-CA630, protease inhibitor, DTT) | IGEPAL works on frozen tissue, Triton-X100 does not (Supp. Fig. 1b) |
| 2 | **High-salt DMF tagmentation buffer** (Tris-acetate, K-acetate, Mg-acetate, 17.6 % DMF), **60 min** | more fragments per nucleus (Supp. Fig. 1c, e) |
| 3 | **EDTA quench before the sort, SDS strip after the sort, Triton-X to neutralise SDS** before PCR | maximal fragment recovery (Supp. Fig. 1d) |
| 4 | **25 nuclei per PCR well**, 11 cycles, calibrated by qPCR | avoid over-amplification (Supp. Fig. 1i) |
| 5 | **8-nt PCR indices** (standard Illumina Nextera N7xx / S5xx + custom X7xx / X5xx) | upstream also notes this |
| 6 | Index reads run **through the constant linkers** (43 + 37 cycles), and a **spike-in library with complemented linkers** balances the bases | the linker is the same in every cluster |

Concept: [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md). There is no RT,
template switching or ligation.

## 2. Oligos

🟢 Verbatim from Supplementary Table 5 (`MOESM5`). No modifications are given except on
pMENTS (`5Phos/` = 5' phosphate).

### 2.1 Tagmentation oligos (sheet "Tagmentation Oligos")

```
p5_1         TCGTCGGCAGCGTCTCCACGCTATAGCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
p5_2         TCGTCGGCAGCGTCTCCACGCATAGAGGCGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
p5_3         TCGTCGGCAGCGTCTCCACGCCCTATCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
p5_4         TCGTCGGCAGCGTCTCCACGCGGCTCTGAGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
p5_5         TCGTCGGCAGCGTCTCCACGCAGGCGAAGGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
p5_6         TCGTCGGCAGCGTCTCCACGCTAATCTTAGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
p5_7         TCGTCGGCAGCGTCTCCACGCCAGGACGTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
p5_8         TCGTCGGCAGCGTCTCCACGCGTACTGACGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
p5_spike-in  AGCAGCCGTCGCAGAGGTGCGTATAGCCTGCGATCGAGGACGGCAGATGTGTATAAGAGACAG
p7_1         GTCTCGTGGGCTCGGCTGTCCCTGTCCCGAGTAATCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
p7_2         GTCTCGTGGGCTCGGCTGTCCCTGTCCTCTCCGGACACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
p7_3         GTCTCGTGGGCTCGGCTGTCCCTGTCCAATGAGCGCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
p7_4         GTCTCGTGGGCTCGGCTGTCCCTGTCCGGAATCTCCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
p7_5         GTCTCGTGGGCTCGGCTGTCCCTGTCCTTCTGAATCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
p7_6         GTCTCGTGGGCTCGGCTGTCCCTGTCCACGAATTCCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
p7_7         GTCTCGTGGGCTCGGCTGTCCCTGTCCAGCTTCAGCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
p7_8         GTCTCGTGGGCTCGGCTGTCCCTGTCCGCGCATTACACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
p7_9         GTCTCGTGGGCTCGGCTGTCCCTGTCCCATAGCCGCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
p7_10        GTCTCGTGGGCTCGGCTGTCCCTGTCCTTCGCGGACACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
p7_11        GTCTCGTGGGCTCGGCTGTCCCTGTCCGCGCGAGACACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
p7_12        GTCTCGTGGGCTCGGCTGTCCCTGTCCCTATCGCTCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
p7_spike-in  CAGAGCACCCGAGCCGACAGGGACAGGCGAGTAATCACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
pMENTS       5Phos/CTGTCTCTTATACACATCT

Read 1 Sequencing Primer    GCGATCGAGGACGGCAGATGTGTATAAGAGACAG
Read 2 Sequencing Primer    CACCGTCTCCGCCTCAGATGTGTATAAGAGACAG
Index 1 Sequencing Primer   CTGTCTCTTATACACATCTGAGGCGGAGACGGTG
```

The Methods call the p5 set "A" and the p7 set "B" transposons (8 × 12 = 96). 🟢

### 2.2 PCR primers (sheet "PCR Primer")

The table lists **49 i7 primers** (N701–N729, X730–X754, X755_spike-in) and **33 i5
primers** (S502–S511, X512–X535, X536_spike-in). Legend in the sheet: "N7.. / S5..
Illumina Nextera index", "X.. Custom index". Representative rows 🟢:

```
N701           CAAGCAGAAGACGGCATACGAGATTCGCCTTAGTCTCGTGGGCTCGG
N702           CAAGCAGAAGACGGCATACGAGATCTAGTACGGTCTCGTGGGCTCGG
X754           CAAGCAGAAGACGGCATACGAGATCGTAATTCGTCTCGTGGGCTCGG
X755_spike-in  CAAGCAGAAGACGGCATACGAGATCTAAGAACCAGAGCACCCGAGCC
S502           AATGATACGGCGACCACCGAGATCTACACCTCTCTATTCGTCGGCAGCGTC
X535           AATGATACGGCGACCACCGAGATCTACACTTACGACCTCGTCGGCAGCGTC
X536_spike-in  AATGATACGGCGACCACCGAGATCTACACTCCTACGGAGCAGCCGTCGCAG
```

The 8-nt indices, cut out of every row with `lib/` 🟡 *(computed)*. For i7 the first
value is as written in the oligo, the second (in brackets) is its reverse complement, which
is what Index 1 reads (and the Illumina "N7xx" name sequence):

```
N701 TCGCCTTA (TAAGGCGA)  N702 CTAGTACG (CGTACTAG)  N703 TTCTGCCT (AGGCAGAA)  N704 GCTCAGGA (TCCTGAGC)
N705 AGGAGTCC (GGACTCCT)  N706 CATGCCTA (TAGGCATG)  N707 GTAGAGAG (CTCTCTAC)  N710 CAGCCTCG (CGAGGCTG)
N711 TGCCTCTT (AAGAGGCA)  N712 TCCTCTAC (GTAGAGGA)  N714 TCATGAGC (GCTCATGA)  N715 CCTGAGAT (ATCTCAGG)
N716 TAGCGAGT (ACTCGCTA)  N718 GTAGCTCC (GGAGCTAC)  N719 TACTACGC (GCGTAGTA)  N721 GCAGCGTA (TACGCTGC)
N722 CTGCGCAT (ATGCGCAG)  N723 GAGCGCTA (TAGCGCTC)  N724 CGCTCAGT (ACTGAGCG)  N726 GTCTTAGG (CCTAAGAC)
N727 ACTGATCG (CGATCAGT)  N728 TAGCTGCA (TGCAGCTA)  N729 GACGTCGA (TCGACGTC)  X730 TACCAGAG (CTCTGGTA)
X731 GGATGGAA (TTCCATCC)  X732 ATTGAGGC (GCCTCAAT)  X733 CGGATAGA (TCTATCCG)  X734 TGGTAGAC (GTCTACCA)
X735 ACCTGGTT (AACCAGGT)  X736 CAGTTCTG (CAGAACTG)  X737 TCGAACGT (ACGTTCGA)  X738 CGTTGCTT (AAGCAACG)
X739 TACCGTTC (GAACGGTA)  X740 TAGGTTGC (GCAACCTA)  X741 GAGGCTAA (TTAGCCTC)  X742 CGACCATA (TATGGTCG)
X743 AGGCAGTA (TACTGCCT)  X744 ATCAAGCG (CGCTTGAT)  X745 CATTGAAG (CTTCAATG)  X746 CGACTTAT (ATAAGTCG)
X747 TCTATACG (CGTATAGA)  X748 AGCATTAG (CTAATGCT)  X749 AATTGGCA (TGCCAATT)  X750 AGATTCGT (ACGAATCT)
X751 TTCATGAC (GTCATGAA)  X752 TGAACTTG (CAAGTTCA)  X753 ATGGCATA (TATGCCAT)  X754 CGTAATTC (GAATTACG)

S502 CTCTCTAT  S503 TATCCTCT  S505 GTAAGGAG  S506 ACTGCATA  S507 AAGGAGTA  S508 CTAAGCCT
S510 CGTCTAAT  S511 TCTCTCCG  X512 TCGACTAG  X513 TTCTAGCT  X514 CCTAGAGT  X515 GCGTAAGA
X516 AAGGCTAT  X517 GAGCCTTA  X518 TTATGCGA  X519 ATCTGAGT  X520 GGATACTA  X521 TAAGATCC
X522 AAGAGATG  X523 AATGACGT  X524 GAAGTATG  X525 ATAGCCTT  X526 TTGGAAGT  X527 ATTCGTTG
X528 AGGATAAC  X529 TTCATCCA  X530 AACGAACG  X531 TGCCTTAC  X532 CGAATTCC  X533 GGTTAGAC
X534 TCCGGTAA  X535 TTACGACC
```

### 2.3 How the oligos interlock: 🟡 (computed against `lib/`)

**Transposons.** All 8 p5 oligos are 63 nt, all 12 p7 oligos 69 nt, and they split into
constant parts and one 8-nt barcode:

```
p5_N (63 nt) = nextera.S5 · TCCACGC · <p5 bc:8> · GCGATCGAGGACGGC · nextera.ME
p7_N (69 nt) = nextera.S7 · CTGTCCCTGTCC · <p7 bc:8> · CACCGTCTCCGCCTC · nextera.ME
```

- **pMENTS = `nextera.ME_RC`** (19 nt, 5'-phosphate): the shared bottom strand of every
  transposon, annealed to the ME at the 3' end of each A/B oligo, the standard
  Tn5-transposon duplex.
- **Read 1 primer = `GCGATCGAGGACGGC` + ME** (34 nt) = the last 34 nt of every p5 oligo;
  **Read 2 primer = `CACCGTCTCCGCCTC` + ME** = the last 34 nt of every p7 oligo. Both
  reads therefore start at the first genomic base, skipping both barcodes.
- **Index 1 primer = revcomp(Read 2 primer)** exactly. No Index 2 primer is listed (§5).
- **PCR primers.** All 48 i7 primers = **`illumina.P7` + 8-nt index + `nextera.S7`**
  (47 nt), all 32 i5 primers = **`illumina.P5` + 8-nt index + `nextera.S5`** (51 nt). So
  the PCR primers prime on the first 15 / 14 nt of the B / A transposon. These are the
  ordinary Nextera index primers (N701 here is the Illumina N701 layout; S502 the S502).
  The i7 index is written as the **reverse complement** of its name (N701 contains
  `TCGCCTTA`, name `TAAGGCGA`); the i5 index is written **as named**.
- **Agreement with the pipeline barcodes.** In `snATAC-master/barcodes/`:
  `r5_ATAC` = the p5 barcodes **as written**, `r7_ATAC` = the p7 barcodes **reverse
  complemented**, `i7_ATAC` = revcomp of the N701–N715 indices (the first 12), `i5_ATAC` =
  S502–S511 as written (the first 8). This is exactly the orientation in which each is
  read (§5). The pipeline lists only one 8 × 12 PCR set; the table has 32 × 48.
- **Spike-ins use the *complement* (not the reverse complement) of every constant
  linker**, and keep the barcodes:
  - p5_spike-in = `complement(S5)` (`AGCAGCCGTCGCAG`) + `complement(TCCACGC)` (`AGGTGCG`)
    + `TATAGCCT` (= p5_1's barcode) + the normal Read 1 site + ME.
  - p7_spike-in = `complement(S7)` (`CAGAGCACCCGAGCC`) + `complement(CTGTCCCTGTCC)`
    (`GACAGGGACAGG`) + `CGAGTAAT` (= p7_1's barcode) + the normal Read 2 site + ME.
  - X536_spike-in = P5 + `TCCTACGG` + `complement(S5)`; X755_spike-in = P7 + `CTAAGAAC`
    + `complement(S7)`: they prime on the spike transposons and not on the regular ones.
  - Effect: in every index cycle that falls on a linker base, the regular library shows
    base X and the spike-in library its complement. A/T and C/G are each one red-laser and one
    green-laser base on the 4-colour HiSeq, so every such cycle has signal in both
    channels. 🟡 (inferred; the paper says only that the spike-in "balances the bases").
    The spike's own index reads stay readable because its read 1/2 sites and ME are the
    normal ones.
- The transposon sequences are **identical** to the Amini 2014 T5/T7 oligos quoted in
  the sci-ATAC-seq note (`TCGTCGGCAGCGTCTCCACGC <8> GCGATCGAGGACGGC ME`, and the p7
  analogue), as the table footnote says. 🟡 (string comparison with that note)

## 3. Step by step

🟢 Online Methods ("Transposome generation", "Single-nucleus ATAC-seq") unless marked.

**Transposomes.** Each A (p5_N) and B (p7_N) oligo is annealed to pMENTS separately
(95 °C 2 min, cooled to 14 °C at 0.1 °C/s), mixed 1:1 molar with unloaded Tn5 (supplied by
Illumina), 30 min at room temperature; then A and B transposomes are combined → 96 A/B
combinations. In a Tn5 loaded with one A and one B, each fragment end gets A or B at
random, so only A…B fragments amplify (§3, step 6). 🟡 (the standard
[Tn5 argument](../ref/concepts/tn5-tagmentation.md))

1. **Nuclei from frozen tissue.** 5–20 mg frozen tissue in 1 mL NPB (5 % BSA, 0.2 %
   IGEPAL-CA630, cOmplete, 1 mM DTT in PBS), 15 min at 4 °C; 30-µm filter; 500 g 5 min.
   (The Reporting summary's pulverisation in liquid N₂ belongs to the separate NeuN-sort
   prep, not to this lysis.) Pellet in
   500 µL 1.1× DMF buffer (36.3 mM Tris-acetate pH 7.8, 72.6 mM K-acetate, 11 mM
   Mg-acetate, 17.6 % DMF); count; 500 nuclei/µL.
2. **Barcoded tagmentation.** 4,500 nuclei per well of a 96-well plate + 1 µL barcoded
   Tn5 transposome (0.25 µM); 60 min at 37 °C, shaking at 500 rpm. Each well = one A/B pair.
3. **Quench.** 10 µL 40 mM EDTA, 37 °C 15 min. Add 20 µL sort buffer (2 % BSA, 2 mM EDTA in
   PBS); pool all wells; 30-µm filter; DRAQ7 3 µM.
4. **Sort.** Sony SH800, **25 nuclei per well** into 4 × 96 = 384 wells containing
   18.5 µL EB with "50 pM Primer i7" and 200 ng BSA. Gating: size, doublet exclusion,
   DRAQ7-high (Supp. Fig. 2).
5. **Strip Tn5.** 2 µL 0.2 % SDS, 55 °C 7 min with shaking; then 2.5 µL 10 % Triton-X to
   quench the SDS (Supp. Fig. 1d).
6. **Indexing PCR.** 2 µL 25 µM Primer i5 + 25 µL NEBNext High-Fidelity 2× master mix;
   **72 °C 5 min** (gap fill-in of the 9-nt Tn5 gap 🟡, not stated in the Methods), 98 °C
   30 s, 11 × (98 °C 10 s, 63 °C 30 s, 72 °C 60 s), hold at 72 °C. The i7 primer's 3' S7
   anneals to the B side, the i5 primer's 3' S5 to the A side; A…A and B…B fragments see
   only one primer and do not amplify exponentially. 🟡
7. **Pool and clean.** All wells (about 15.5 mL) + 80 mL PB with pH indicator + 4 mL 3 M
   sodium acetate pH 5.2; four MinElute columns, 15 µL EB each, pooled. Double-sided
   AMPure XP: 0.55× (keep the supernatant), then 1.5× more beads; 2 × 80 % EtOH; elute in
   20 µL EB. Qubit, Tapestation (nucleosomal ladder).

**Species-mixing variant** (E15.5 forebrain + GM12878, 1:1 before tagmentation): only 96
wells after sorting, **13** PCR cycles, MiSeq (§5). 🟢

**Bulk ATAC-seq** (Methods "ATAC-seq", after "NeuN-negative sorting": 20,000 sorted
nuclei, 0.5 µM Tn5 in DMF buffer, 60 min, MinElute, 8–10 PCR cycles) 🟢 is a separate
assay and not part of the single-nucleus library.

## 4. Final library: 🟡 (assembled from the oligos above)

Top strand, 5'→3', the P5 end first; `<…>` = variable. Lengths are from `lib/`.

```
P5 (29) · <i5:8> · S5 (14) · TCCACGC (7) · <p5 bc:8> · GCGATCGAGGACGGC (15) · ME (19)
  · <genomic insert>
  · ME_RC (19) · GAGGCGGAGACGGTG (15) · <revcomp p7 bc:8> · GGACAGGGACAG (12) · S7_RC (15)
  · <revcomp i7-as-written:8> · P7_RC (24)
```

Constant length outside the insert: 29 + 8 + 14 + 7 + 8 + 15 + 19 + 19 + 15 + 8 + 12 + 15 +
8 + 24 = **201 nt**. 🟡 *(computed)*

The cell barcode is the 4-tuple (p5 bc, p7 bc, i5, i7): 96 tagmentation combinations ×
the PCR combinations used. The table has 32 i5 × 48 i7; Supp. Table 1 says the two
replicates were sequenced together and told apart by "replicate specific barcode
combinations (Set_1 or Set_2)" 🟢. So each replicate probably used half the i5 or i7
set, e.g. 16 × 24 = 384 wells. 🟡 (inferred; the split is not stated 🔴)

**Spike-in library** (same layout): `complement(S5) · AGGTGCG` replaces `S5 · TCCACGC`,
`GACAGGGACAGG · complement(S7)` replaces `CTGTCCCTGTCC · S7` (on the top strand,
reverse-complemented accordingly), and the PCR indices are `TCCTACGG` (i5) /
`CTAAGAAC` as written (i7). 🟡

## 5. Read layout / sequencing

🟢 HiSeq 2500, 25 pM per lane, **custom primers** (Read 1, Read 2, Index 1 from
Supplementary Table 5), **50 + 43 + 37 + 50** (Read 1 + Index 1 + Index 2 + Read 2).
"The first 8 bp of Index1 correspond to the p7 barcode and the last 8 bp to the i7
barcode. The first 8 bp of Index2 correspond to the i5 barcode and the last 8 bp to the
p5 barcode." Mixing experiment: MiSeq, 15 pM, 44 + 43 + 37 + 44.

Computed from the structure 🟡 *(computed)*, in agreement with the stated lengths:

| Read | Primer | Reads | Length |
|---|---|---|---|
| Read 1 | custom, `GCGATCGAGGACGGC`+ME | genomic insert from the A end | 50 |
| Index 1 | custom = revcomp(Read 2 primer) | `<revcomp p7 bc:8>` · `GGACAGGGACAG` · `S7_RC` · `<i7 name:8>` = 8 + 27 + 8 | **43** ✓ |
| Index 2 | none listed → standard, from the P5 end | `<i5:8>` · `S5` · `TCCACGC` · `<p5 bc:8>` = 8 + 21 + 8 | **37** ✓ |
| Read 2 | custom, `CACCGTCTCCGCCTC`+ME | genomic insert from the B end | 50 |

- The 27 + 21 linker bases are **read in light cycles**, not skipped with dark cycles as
  in Amini 2014 / sci-ATAC-seq 2018. Hence the spike-in. 🟡
- Index 2 reads **i5 first, then the Tn5 barcode**: the HiSeq 2500 forward-strand order,
  the same as Amini 2014 on the HiSeq 2000. On a NextSeq / reverse-complement workflow
  with a custom Index 2 primer (as sci-ATAC-seq 2018) the order flips. 🟡
- The pipeline's barcode files (`r5`, `r7`, `i5`, `i7`) are in exactly the orientation
  this layout predicts (§2.3), which confirms it independently of the text. 🟡 *(computed)*
- 8-nt barcodes; the Supp. Table 1 legend notes that for E11.5 "7 out of 8 bp were
  detected for the p5 barcode", i.e. Index 2 was one cycle short in that run. 🟢

## 6. Check against upstream (scg_lib_structs)

The upstream page is short; it refers to its sci-ATAC-seq page for the steps.

| Upstream says | Paper / supplement | Verdict |
|---|---|---|
| same combinatorial strategy and oligo design as sci-ATAC-seq | transposon + sequencing-primer sequences are the Amini 2014 ones (table footnote) | agree 🟢 |
| modifications for frozen tissue | NPB lysis, DMF buffer, EDTA → SDS → Triton (Supp. Fig. 1) | agree 🟢 |
| slightly different PCR primers, 8-bp index | P5/P7 + 8-nt index + S5/S7 | agree 🟢 |
| spike-in oligos to balance bases in the middle of the index reads | yes; computed to be complemented linkers | agree, mechanism added 🟡 |
| "step-by-step … same as sci-ATAC-seq" | differs in the sort (25 nuclei/well into 384 wells), SDS/Triton, single 11-cycle PCR, HiSeq 2500 read order (i5 before p5 in Index 2), light linker cycles | partly; upstream gives no read lengths for snATAC-seq |

## 7. Open questions

- 🔴 **"50 pM Primer i7"** in the sort buffer, against 2 µL of 25 µM i5 (about 1 µM final).
  50 pM i7 would be far too little for 11 cycles; probably a unit typo (nM or µM?). Not
  resolvable from the sources.
- 🔴 Which i5 × i7 combinations make up "Set_1" and "Set_2", and why 32 i5 and 48 i7 are
  listed for 384 wells.
- 🔴 How much spike-in library was mixed in, and at which step: tagmented and amplified
  separately, or as nuclei? The Methods give no amount.
- 🟡 The spike-in's barcodes equal p5_1 / p7_1, and its PCR indices (`TCCTACGG`,
  `CTAAGAAC`) are not in the regular sets. So the demultiplexer can tell spike reads by
  the PCR index, but this is not stated.
- 🟡 Barcoded Tn5 tagments less efficiently and gives larger fragments (about 550 vs 300 bp,
  Supp. Fig. 1g); the 0.25 µM loading is the compromise from Supp. Fig. 1h.

## 8. How this note was made (tool evaluation)

`tools/get_sources.py` fetched the PMC full text, all six Springer supplements (the zip
was unpacked to `moesm6_unzipped/`) and the parent sci-ATAC-seq paper. The four PMC
NIHMS supplement copies were behind the download gate. `tools/doctext.py` turned the
xlsx into tab-separated text that was easy to parse. `tools/scrape_primers.py` put
Supplementary Table 5 first and tagged P7 / S7 / P5 / S5 / ME, but the output repeats a
long `context` window for every row of a dense table. The methods were read by hand. The
barcode orientation was settled with the pipeline's barcode files, which no tool looks at.
