# PIP-seq — droplet scRNA-seq by vortexing barcoded hydrogel beads

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. Relationships marked 🟡 *(computed)* were worked out with `lib/` while
> writing this note; they are not yet asserted in a self-test, because this protocol has
> no `tools/` module yet (status `notes`). Claims taken only from the upstream
> scg_lib_structs page are 🟡 *(upstream)*: a careful secondary source that reconstructed
> parts of the design from public FASTQ files, not a primary one.

**PIP-seq** (particle-templated instant partition sequencing) — Clark IC, Fontanez KM,
Meltzer RH, Xue Y, Hayford C, May-Zhang A, D'Amato C, Osman A, Zhang JQ, Hettige P,
Ishibashi JSA, Delley CL, Weisgerber DW, Replogle JM, Jost M, Phong KT, Kennedy VE,
Peretz CAC, Kim EA, Song S, Karlon W, Weissman JS, Smith CC, Gartner ZJ, Abate AR.
"Microfluidics-free single-cell genomics with templated emulsification."
*Nature Biotechnology* 41, 1557 (2023).
doi:[10.1038/s41587-023-01685-z](https://doi.org/10.1038/s41587-023-01685-z)
(PMC10635830, CC-BY 4.0). Preprint: bioRxiv
doi:[10.1101/2022.06.10.495582](https://doi.org/10.1101/2022.06.10.495582) (v2, 25 Oct 2022).
Data: GEO SuperSeries GSE202919.

Other papers in the catalogue for this protocol:

- Delley CL, Abate AR. "Modular barcode beads for microfluidic single cell genomics."
  *Sci. Rep.* 11, 10857 (2021). doi:[10.1038/s41598-021-90255-x](https://doi.org/10.1038/s41598-021-90255-x)
  — the ligation-based split-pool bead chemistry PIP-seq says it follows (ref. 30 of the
  paper). Its handles are **different** from PIP-seq's (see §1).
- Plessy C. et al. "Linking promoters to functional transcripts in small samples with
  nanoCAGE and CAGEscan." *Nat. Methods* 2010,
  doi:[10.1038/nmeth.1470](https://doi.org/10.1038/nmeth.1470) — cited by upstream only
  for "single-primer semi-suppressive PCR"; not fetched or read here.

Sources read (fetched by `tools/get_sources.py`, into
`$CHEM_DATA/sources/pip-seq__10.1038+s41587-023-01685-z/`, never committed):

| File | What | Used for |
|---|---|---|
| `s41587-023-01685-z_PMC10635830.html.txt` | journal paper, full text (PMC) | **methods** (bead fabrication, RT, WTA, library), read lengths |
| `s41587-023-01685-z_41587_2023_1685_MOESM2_ESM.xlsx.txt` | Supplementary Table 1 (= preprint `media-1.xlsx`, same md5) | **every published oligo** |
| `s41587-023-01685-z_41587_2023_1685_MOESM7_ESM.xlsx.txt` | Supplementary Table 6 (= preprint `media-6.xlsx`) | read lengths per run |
| `2022.06.10.495582_v2.full.pdf.txt` | the preprint | cross-check of the methods (identical in the parts used) |
| `upstream_PIP-seq.html.txt` | scg_lib_structs PIP-seq page (V2 and Fluent V3/V4) | second source, checked below |
| `upstream_PIP-seq_v1p.html.txt` | scg_lib_structs "PIP-seq V1 prototype" page | the obsolete first library type (fetched by hand, see §8) |
| `upstream_data_*.tsv` | upstream barcode whitelists (V2 rounds 1–3, Fluent V3 rounds 1–4) | barcode lengths only (fetched by hand) |
| `s41598-021-90255-x_PMC8149635.html.txt`, `…_MOESM4_ESM.xlsx.txt` | Delley 2021 and its oligo tables | what the bead concept builds on |
| `…MOESM1_ESM.pdf` | Nature "reporting summary" (image-only PDF; read as page images) | nothing about chemistry |
| `…MOESM3/4/5/6/8`, `media-2..5,7` | marker-gene tables, patient table | not chemistry |

Not obtainable / not published:

| What | Status |
|---|---|
| Sequences of the PIP-seq split-pool **barcode and splint oligos** (the ligation blocks between the acrydite primer and the UMI) | 🔴 not in Supplementary Table 1; upstream reconstructed them from reads |
| Fluent Biosciences kit oligos (TSO/WTA primer/adapter/index primers as shipped, user guides FB0002079 etc.) | 🔴 proprietary; upstream gives guesses |
| The PMC copies of the supplements (`(manual)` rows in `MANIFEST.tsv`, JavaScript download gate) | not needed: the same files were fetched from Springer |

---

## 1. What it is

A **droplet, poly(T)-bead, 3'-end** scRNA-seq method (the Drop-seq / inDrop / 10x family)
whose novelty is entirely physical: droplets are made **without a microfluidic chip**. 🟢

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | **Particle-templated emulsification**: barcoded polyacrylamide beads, cells, proteinase K and DTT are mixed, oil is added and the tube is **vortexed**; each droplet forms around one bead, so droplet size is set by bead size and almost every droplet holds a bead 🟢 | Hatori 2018 (ref. 29); HyDrop / inDrop use chips for the same beads |
| 2 | **Heat-triggered lysis**: proteinase K is nearly inactive at 4 °C, so cells survive mixing; 65 °C after emulsification lyses them inside their droplet 🟢 | — |
| 3 | **RT after breaking the emulsion**, in bulk, with mRNA already hybridised to the bead poly(T) 🟢 | inDrop / 10x do RT inside droplets; Drop-seq does it in bulk after breaking, as here |
| 4 | Beads built by **splinted-ligation split-pool** (four rounds of 96) rather than phosphoramidite split-pool 🟢 | Delley & Abate 2021; SPLiT-seq-style ligation |

Downstream of the beads the chemistry is standard: oligo-dT priming,
[template switching](../ref/concepts/template-switching.md) with a SMART-handle TSO,
single-primer whole-transcriptome amplification (WTA), then a 3'-end library made either by
**Nextera XT** ([Tn5](../ref/concepts/tn5-tagmentation.md)) with a custom P5 primer or by
**adapter ligation** (Watchmaker kit). See also
[reverse transcription](../ref/concepts/reverse-transcription.md).

**Relation to Delley 2021** 🟡 *(computed)*: Delley's bead primer (`pBB1`, Supplementary
Data Table 7 of that paper) starts `/5Acryd/ ACTAACAATAAGCTC UAU …` with a T7 promoter and a
custom "PE2c" read-1 handle; none of that is in the PIP-seq acrydite primer, which carries
the SMART handle and a partial TruSeq Read 1 instead. PIP-seq keeps the *method* (acrylamide
beads, T4-ligase split-pool with short splints) and changes the *sequences* — as upstream
also notes.

**Library generations.** The paper describes one method but its data contain (at least)
three read structures, and the commercial kit a fourth. Only one of them is fully backed by
the paper's methods:

| Version | Barcode | Library prep | Primary evidence | Upstream |
|---|---|---|---|---|
| V1 prototype | 2 × 8-nt barcode + inDrop W1 linker, 6-nt UMI | Nextera XT | 🔴 none (only in FASTQ, e.g. the gefitinib runs) | `PIP-seq_v1p` page |
| **V2 (paper methods)** | 3 × 8 nt with 7-nt linkers, 12-nt UMI, T19V | Nextera XT + `PIPs_P5library` + N7xx | 🟢 acrydite primer, UMI/poly(T), TSO, WTA, P5 primer; 🔴 barcode blocks | main page, "PIP-seq V2" |
| Fluent PIPseq V3.0 | 8 + 6 + 6 + 8 nt, 12-nt UMI, T30V | adapter ligation, dual-indexed | 🟡 the paper only says the library was made by adapter ligation with a Watchmaker Genomics kit, and that beads carry "four barcodes" | main page, "FluentBio" |
| Fluent PIPseq V4.0 | as V3 + 0–3 nt phase spacer before barcode 1 | as V3 | 🔴 not in the paper | main page |

## 2. Oligos

🟢 Verbatim from Supplementary Table 1 (`MOESM2_ESM.xlsx`, sheet "Oligonucleotides"; the
preprint's `media-1.xlsx` is byte-identical). Notation as given: `rG` = ribo-G,
`*` = phosphorothioate, `[N5xx]` / `[N7xx]` = index placeholders.

```
PIPS_TSO          AAGCAGTGGTATCAACGCAGAGTGAATrGrGrG
PIPS_WTA_primer   AAGCAGTGGTATCAACGCAGAGT
PIPs_P5library    AATGATACGGCGACCACCGAGATCTACACGCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGT*A*C
P5-PE1            AATGATACGGCGACCACCGAGATCTACAC[N5xx]ACACTCTTTCCCTACACGAC*G*C
Weissman_U6       CAAGCAGAAGACGGCATACGAGAT[N7xx]GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGTGTTTTGAGACTATAAGTATCCCTTGGAGAACCACCTTGT*T*G
```

🟢 In the methods text (not in the table), "Barcode bead fabrication":

```
acrydited primer   /5Acryd/TTTTTTTAAGCAGTGGTATCAACGCAGAGTACGACTCCTCTTTCCCTACACGACGCTCTTCC
UMI with poly(T)   NNNNNNNNNNNNTTTTTTTTTTTTTTTTTTTV
```

The table also lists 19 TotalSeq-A ADT barcodes and two hashtag (HTO) barcodes (15 nt
each, BioLegend) used in the MPAL and hashing experiments — antibody barcodes, not
protocol oligos.

Standard kit primers used but not listed: Nextera XT kit, "standard Nextera P7 indexing
primers (N70x)" 🟢 (sequence not given in the paper; upstream writes the N7xx as
`CAAGCAGAAGACGGCATACGAGAT[8-bp i7 index]GTCTCGTGGGCTCGG`, i.e. P7 + i7 + s7 🟡).

### Upstream's additional oligos — 🟡 (upstream, not in the paper)

For V2 upstream adds a bead splint `pBB2` and three ligation blocks, reconstructed from
Read 1 of SRR19180490:

```
pBB2         /5Phos/AGATCGGAAGAGCGTCGTGTAGGGAAAGAGGAGTCGTACTCTGCGTTGATACCACTGCTT
plate-1-BC   /5Phos/GATCT[8-bp barcode1]ATGCATC          splint plate-1-SP /5Phos/[barcode1 rc]
plate-2-BC   /5Phos/CTCGAGG[8-bp barcode2 rc]GATGCAT     splint plate-2-SP /5Phos/[barcode2]
plate-3-BC   /5Phos/CCTCGAG[8-bp barcode3][12-bp UMI]TTTTTTTTTTTTTTTTTTTV   splint plate-3-Sp /5Phos/[barcode3 rc]
```

For Fluent V3/V4 upstream gives (from Fluent part numbers and example data, explicitly
"educated guesses"): a two-primer WTA (`Forward CTCTTTCCCTACACGACGCTC`, `Reverse` =
PIPS_WTA_primer), a ligation "Library Adapter Mix" (`/5Phos/CTGTCTCTTATACACATCTCCGAGCC`
annealed to a short `TGACAGAGAATAT`), and dual-index PCR primers
`AATGATACGGCGACCACCGAGATCTACAC[i5]ACACTCTTTCCCTACACGACGC` and
`CAAGCAGAAGACGGCATACGAGAT[i7]GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG`.

### How they interlock — 🟡 (computed against `lib/`)

One 23-nt handle, **`rt.SMART_HANDLE`** (`AAGCAGTGGTATCAACGCAGAGT`, the SMART-seq ISPCR
sequence), sits at both ends of the cDNA:

- **Acrydite bead primer** (62 nt) = T₇ spacer + **SMART handle** + `ACGACTC` +
  `CTCTTTCCCTACACGACGCTCTTCC` (= `TRUSEQ_READ1[3:28]`). After the first ligation block
  adds `GATCT` (upstream's plate-1-BC; 🔴 in the paper) the bead carries
  `TRUSEQ_READ1[3:]`, the last **30 of 33 nt** of the TruSeq Read 1 primer site; the first
  three (`ACA`) are replaced by `CTC`. 🟡
- **PIPS_TSO** (30 nt) = **SMART handle** + `GAAT` + `rGrGrG`. 🟢 sequence; the `GAAT`
  spacer is not explained 🔴. After template switching every first strand ends
  `…CCC ATTCACTCTGCGTTGATACCACTGCTT-3'` (revcomp of the TSO DNA part).
- **PIPS_WTA_primer** is exactly the **SMART handle** (23 nt), so one primer amplifies
  from both ends — single-primer PCR, as in SMART-seq2. Short inserts form panhandles
  (both ends complementary) and are suppressed. 🟡
- **PIPs_P5library** (66 nt) = **`illumina.P5`** + 12-nt `GCCTGTCCGCGG` (unexplained 🔴) +
  **SMART handle** + `AC`, with the last two bases phosphorothioate-protected. The `AC`
  matches the bead side (handle + `ACGACTC…`) but not the TSO side (handle + `GAAT…`), so
  this primer extends **only from the barcode end** of the cDNA. That is how a library made
  with a shared handle at both ends still becomes a 3'-end library. 🟡 (computed; upstream
  draws the same conclusion)
- **P5-PE1** (used only for the CROP-seq sgRNA library) = **`illumina.P5`** + i5 +
  `TRUSEQ_READ1[:22]` (`ACACTCTTTCCCTACACGACGC`). Its 3' end overlaps the bead's
  `TRUSEQ_READ1[3:]` by 19 nt and restores the full Read 1 site. The same layout is
  upstream's Fluent "Library P5 Index", so the paper's table is primary evidence for that
  part of the V3 design. 🟡
- **Weissman_U6** = **`illumina.P7`** + i7 + **`nextera.ADAPTOR_S7`** (s7 + ME, 34 nt) +
  a 43-nt sgRNA-cassette-specific 3' part (`GTGTTTTGAG…CCTTGT*T*G`). So the sgRNA amplicon
  is read with the **Nextera Read 2 / i7 primers**, the same as the Nextera XT transcriptome
  library, and both can share a run. Where the 43 nt bind in the CROP-seq vector was not
  checked 🔴.

Upstream's oligos 🟡 *(computed)*:

- `pBB2` is the exact reverse complement of acrydite primer (without T₇) + `GATCT`: a splint
  that pairs with the whole handle–Read 1 region and leaves a 5' `AGATC` overhang for
  the first block. Plausible, but the paper only says the beads were annealed to "a
  complementary oligonucleotide".
- V3 WTA Forward = `TRUSEQ_READ1[3:24]`; the V3 P5 index tail = `TRUSEQ_READ1[:22]`.
- V3 "Library Adapter Mix" top strand = `nextera.ME_RC` + `S7_RC[:7]` (`CTGTCTCTTATACACATCT`
  + `CCGAGCC`); the P7 index primer tail is exactly `nextera.READ2_PRIMER`. So even the
  ligation library is read with the **Nextera** Read 2 and Index 1 primers. 🟡

## 3. Bead fabrication (before the experiment)

Conditions 🟢 from "Barcode bead fabrication" (prototype beads; later commercial beads
from Fluent BioSciences, e.g. FB0003067 / FB0002617).

1. **Beads**: co-flow drop maker, 6 % acrylamide / 0.1 % bis, APS, TEMED in the oil, 50 µM
   acrydite primer; polymerised 12 h at room temperature; 80 µm beads. The primer is
   copolymerised through its 5' acrydite.
2. **Split-pool ligation**, per round: beads in T4 ligase buffer, heated with a
   complementary oligo to 75 °C 2 min and cooled (anneal); 100 µL beads per well of a
   96-well plate holding a unique barcode, 1.9 U/µL T4 DNA ligase, 25 °C 1 h, 65 °C 10 min;
   pool, wash 5 × TET. Repeated "to add four barcodes and a UMI with poly(T)". 🟢
3. The main text says the beads were made by split-pool ligation of four 6-bp randomers, giving ~10⁸ (96⁴) barcodes. 🟢
   This does **not** match V2 (three 8-nt barcodes = 96³) and only partly matches V3
   (8 + 6 + 6 + 8 nt); see §7.
4. Make single-stranded (upstream: NaOH strip of the splints 🟡; in Delley 2021, 100 mM NaOH
   washes 🟢 for that paper).

Bead oligo after assembly, V2 🟡 (upstream; constant parts 🟢 where marked):

```
5'- /5Acryd/ · T7 · SMART handle (🟢) · ACGACTC · TruSeq R1[3:28] (🟢) · GATCT · <bc1 8> · ATGCATC · <bc2 8> · CCTCGAG · <bc3 8> · <UMI 12> · T19 · V -3'
```

The `<UMI 12> · T19 · V` tail is 🟢 (methods). Fluent V3 🟡 (upstream):

```
5'- /5Acryd/ · T7 · SMART handle · ACGACTC · TruSeq R1[3:28] · GATCT · <bc1 8> · ATG · <bc2 6> · GAG · <bc3 6> · TCGAG · <bc4 8> · <UMI 12> · T30 · V -3'
```

## 4. Step by step (V2, the paper's methods)

Volumes and conditions 🟢 from "PIP-seq workflow" (standard small-format run).

1. **Mix**: 5 µL cells (500/µL, in PBS + 1 % Pluronic F127) into 35 µL barcoded beads
   carrying 29 U/mL proteinase K and 70 mM DTT; 10 pipette strokes.
2. **Emulsify**: 280 µL fluorinated oil with surfactant; vortex 15 s horizontal + 2 min
   vertical at maximum speed. Remove 230 µL excess oil.
3. **Lyse**: 65 °C 35 min, 4 °C hold. mRNA poly(A) anneals to the bead poly(T) in the
   droplet. The emulsion can be held at 0 °C for 72 h at this point without loss of quality.
4. **Break**: 180 µL high-salt buffer (250 mM Tris pH 8, 375 mM KCl, 15 mM MgCl₂, 50 mM DTT)
   + 40 µL perfluoro-octanol; wash beads 3 × in 2× RT buffer with 1 % Pluronic F68. The salt
   keeps mRNA hybridised to the beads once the droplets merge. 🟡
5. **RT with template switching** (bulk, on beads): 25 µL beads + 25 µL mix — 4.8 % PEG8000,
   4 % PM400, 2.5 µM PIPS_TSO, 1 mM dNTPs, RNase inhibitor, Maxima H-minus RT (1 U/µL);
   25 °C 30 min, 42 °C 90 min, 85 °C 10 min. First strand, still on the bead 🟡:
   `5'-bead oligo · (dT)V · cDNA · CCC · revcomp(TSO DNA)-3'` — see
   [template switching](../ref/concepts/template-switching.md).
6. **WTA** without purification: + 50 µL 2× KAPA HiFi, 0.25 µM PIPS_WTA_primer; 95 °C 3 min;
   16 × (98 °C 15 s, 67 °C 20 s, 68 °C 4 min); 72 °C 5 min. The first cycle copies the
   bead-bound strand into solution; thereafter both ends carry the SMART handle. 🟡
7. **Clean-up**: beads removed on a Spin-X filter; 0.6× AMPure XP.
8. **Library (Nextera XT)**: tagment the WTA product; amplify with **PIPs_P5library** and
   a standard **N70x** primer. Only fragments with the barcode end on one side and Tn5's
   s7 adapter on the other amplify: PIPs_P5library primes only on the bead end (§2) and the
   N7xx only on s7. 🟡 (cycle numbers and input amount not given 🔴)
9. **Sequence**: NextSeq 2000, 15 % PhiX.

Variations in the same paper 🟢:

- **Adapter-ligation libraries** (72-h hold, breast tissue, high-cell-count PBMC): "performed according to manufacturer's instructions (Watchmaker Genomics,
  7K0019-024)". No further detail; this is presumably the Fluent V3 route upstream
  reconstructs. 🟡
- **Large formats** (e.g. the high-cell-number breast study): 800 µL beads + 40 µL cells + 4 mL oil; lysis then used
  a separate "lysis emulsion" (FB0003039) at 37 °C 45 min instead of 65 °C.
- **CROP-seq sgRNA library**: 1 ng WTA cDNA, 0.5 µM P5-PE1 + 0.25 µM Weissman_U6, KAPA HiFi;
  95 °C 3 min; 10 × (95 °C 20 s, 70 °C 30 s −0.2 °C/cycle, 72 °C 20 s); 8 × (95 °C 20 s,
  68 °C 30 s, 72 °C 20 s); 72 °C 4 min; 0.5×/0.8× AMPure. Pooled 20:1 with the transcriptome
  library.
- **ADT / HTO** (TotalSeq-A): "processed according to the PIP-seq Single Cell Epitope
  Sequencing user guide (FB0002079)". The ADT library chemistry is not described. 🔴

## 5. Final libraries — 🟡 (assembled from the oligos above)

**V2, Nextera XT** (the paper's methods; constant parts computed, barcode linkers from
upstream):

```
5'- P5 · GCCTGTCCGCGG · SMART handle · ACGACTC · TruSeq R1[3:] · <bc1 8> · ATGCATC · <bc2 8> · CCTCGAG · <bc3 8> · <UMI 12> · (dT)V · <cDNA> · ME' · s7' · <i7 8> · P7' -3'
```

The Read 1 primer anneals to the 30-nt `TruSeq R1[3:]`; its first three bases (`ACA`) sit
opposite `CTC` and do not pair. That is unusual but enough for priming. 🟡 *(computed)*.
Upstream's drawing leaves 4 nt unpaired; 3 is the computed number.

**Fluent V3, adapter ligation** 🟡 (upstream):

```
5'- P5 · <i5 8> · TruSeq R1 · <bc1 8> · ATG · <bc2 6> · GAG · <bc3 6> · TCGAG · <bc4 8> · <UMI 12> · (dT)V · <cDNA> · A · ME' · s7' · <i7 8> · P7' -3'
```

Here the Read 1 site is complete (restored by the P5 index primer tail). V4 inserts 0, 1, 2
or 3 nt (`—`, `T`, `GT`, `TGA`) between TruSeq R1 and bc1, so that clusters do not all read
the same base in the same cycle. 🟡 (upstream, from Fluent example FASTQ)

**V1 prototype** 🟡 (upstream page only): bead oligo
`… TruSeq R1[3:] · GATCT · <0–3 nt spacer> · <bc1 8> · GAGTGATTGCTTGTGACGCCTT (inDrop W1) · <bc2 8> · <UMI 6> · T19V`,
amplified with a **different** "PIPs_P5library",
`AATGATACGGCGACCACCGAGATCTACACTAGATCGCCTCTTTCCCTACACGACGC` = P5 + `TAGATCGC` +
`TRUSEQ_READ1[3:22]` *(computed)*. That sequence is not in the paper's Supplementary
Table 1. It is upstream's reading of the data, and the paper may have reused one name for
two primers. 🔴

## 6. Sequencing

🟢 Illumina NextSeq 2000. Read lengths from Supplementary Table 6 (Read 1 / Read 2):
51/70 (2022 runs: species mix, PBMC, breast T200, hashing), 51/66 (breast re-sequencing),
55/66 (72-h hold, MPAL, 2021 species mix), 55/91 (gefitinib), 55/145 (CROP-seq).

🟡 *(computed)* Read 1 needed for barcodes + UMI: V2 = 8+7+8+7+8+12 = **50 nt**; V3 =
8+3+6+3+6+5+8+12 = **51 nt** — exactly the 51-cycle Read 1 of the 2022 runs. V1 needs up
to 3+8+22+8+6 = 47.

Primers 🟡 (upstream; consistent with the oligos above): Read 1 = `illumina.TRUSEQ_READ1`;
Index 1 = `nextera.INDEX1_PRIMER`; Read 2 = `nextera.READ2_PRIMER` (cDNA, starting at the
fragment end far from poly(A)). V3 adds an i5 read with
`illumina.INDEX2_PRIMER_RC` (revcomp of TruSeq Read 1). V2 has no i5 index:
PIPs_P5library has none.

## 7. Open questions and disagreements

- 🔴 **Barcode geometry.** The paper says four ligated barcodes and "four 6-base pair (bp)
  randomers" (~96⁴). Upstream's V2 whitelists are three rounds of 96 × 8 nt; its V3 whitelists are four
  rounds of 8/6/6/8 nt (checked: `upstream_data_*.tsv`). Neither is "four 6-bp". The
  methods' poly(T) (`T19V`) agrees with upstream V2, **not** with upstream V3 (`T30V`).
  The paper therefore seems to describe V2 beads with a V3-like barcode count.
- 🔴 The **ligation block sequences** (linkers `ATGCATC`, `CCTCGAG`, `ATG`, `GAG`, `TCGAG`,
  the `GATCT` start) and splints are unpublished; upstream inferred them from reads.
- 🔴 Function of `GCCTGTCCGCGG` in PIPs_P5library and of `GAAT` in PIPS_TSO.
- 🔴 Whether the Watchmaker ligation libraries use upstream's Fluent adapter and index
  primers, and whether WTA for those runs is single-primer (paper) or two-primer
  (upstream V3).
- 🔴 Which "PIPs_P5library" made the V1-prototype (gefitinib) libraries (§5).
- 🔴 ADT/HTO library chemistry (Fluent user guide FB0002079 not public).
- 🟡 The 43-nt 3' part of Weissman_U6 — not matched to a vector sequence here.

## 8. How this note was made (tool evaluation)

`tools/get_sources.py` fetched the PMC full text, all eight Springer supplements, the
bioRxiv v2 preprint with supplements, and Delley 2021 with its supplements; the PMC
supplement copies were behind the download gate but were not needed. The upstream page
links two further pages it does not fetch: the **V1-prototype page**
(`methods_html/PIP-seq_v1p.html`) and the **barcode whitelists**; both were fetched by hand
into the sources directory. `MOESM1_ESM.pdf` (Microsoft "Print to PDF", image-only) gave an
empty text twin; rendered to images, it is only the reporting summary. `scrape_primers.py`
found every oligo in Supplementary Table 1 and the acrydite primer in the PMC text with its
`/5Acryd/`. In the preprint, a page line number breaks the acrydite primer and the scraper
returned it 10 nt short. Checking (`--find`): oligos whose 3' end carries
phosphorothioate stars (`…GT*A*C`, `…CGAC*G*C`, `…TTGT*T*G`) are only found in the
supplementary table when the starred tail is left off; the full stripped sequence matched
only upstream's copy (or nothing, for Weissman_U6). The table rows were checked by eye.
The 15-nt ADT/HTO barcodes and the Delley barcode tables made up
most of the hits.
