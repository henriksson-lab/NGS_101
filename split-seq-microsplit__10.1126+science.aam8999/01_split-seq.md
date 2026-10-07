# SPLiT-seq and microSPLiT — single-cell RNA-seq by split-pool ligation barcoding

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published, or published only in a file that could not be fetched. Relationships marked
> 🟡 *(computed)* were worked out with `lib/` while writing this note; they are not yet
> asserted in a self-test, because this protocol has no `tools/` module yet (status
> `notes`). Claims taken only from the upstream scg_lib_structs page are 🟡 *(upstream)*:
> a careful secondary source, not the paper.

**SPLiT-seq** — Rosenberg AB, Roco CM, Muscat RA, Kuchina A, Sample P, Yao Z, Gray L,
Peeler DJ, Mukherjee S, Chen W, Pun SH, Sellers DL, Tasic B, Seelig G. "Single-cell
profiling of the developing mouse brain and spinal cord with split-pool barcoding."
*Science* 360, 176–182 (2018). doi:[10.1126/science.aam8999](https://doi.org/10.1126/science.aam8999)
(PMC7643870, author manuscript). Preprint: "Scaling single cell transcriptomics through
split pool barcoding", *bioRxiv* 2017, doi:[10.1101/105163](https://doi.org/10.1101/105163)
— an earlier chemistry (see §7).

**microSPLiT** — Kuchina A, Brettner LM, Paleologu L, Roco CM, Rosenberg AB, Carignano A,
Kibler R, Hirano M, DePaolo RW, Seelig G. "Microbial single-cell RNA sequencing by
split-pool barcoding." *Science* 371, eaba5257 (2021).
doi:[10.1126/science.aba5257](https://doi.org/10.1126/science.aba5257) (PMC8269303).
Preprint: *bioRxiv* 2019, doi:[10.1101/869248](https://doi.org/10.1101/869248) (v2, 11 Dec
2019, CC-BY-NC-ND 4.0). A later step-by-step protocol, *Nature Protocols* 2024,
doi:10.1038/s41596-024-01007-w, is catalogued separately and was **not** read for this note.

Sources read (in `$CHEM_DATA/sources/split-seq-microsplit__10.1126+science.aam8999/`,
never committed):

| File | What | Used for |
|---|---|---|
| `upstream_SPLiT-seq.html.txt` | scg_lib_structs page (SPLiT-seq + microSPLiT) | the only full account of the **published SPLiT-seq** oligos, steps, library, read layout |
| `science.aam8999_PMC7643870.html.txt` (+ `_efetch.xml.txt`) | SPLiT-seq main text (PMC) | concept, barcode arithmetic; **no methods** (they are in the gated supplement) |
| `105163_v1.full.pdf.txt` | SPLiT-seq preprint v1 | the earlier, RT-then-three-ligations design (§7) |
| `869248_v2.full.pdf.txt` | microSPLiT preprint v2 | concept, sequencing summary |
| `869248_media-1.pdf.txt` | microSPLiT preprint Supplementary Materials | **all microSPLiT methods**, incl. library prep |
| `869248_media-2.xlsx.txt` | microSPLiT preprint Auxiliary Table S3 | **every microSPLiT oligo**, all 96 × 3 plate barcodes |
| `science.aba5257_PMC8269303_efetch.xml.txt` | published microSPLiT text (NCBI efetch; the PMC HTML was a captcha page) | "Materials and Methods Summary" — same RT/ligation/lysis conditions, **no library prep** |
| `869248_media-3.csv` | cluster marker genes | not chemistry |

Could not be fetched (PMC download gate / reCAPTCHA; get by hand):

| What | URL |
|---|---|
| SPLiT-seq Supplementary Materials (all methods) | https://pmc.ncbi.nlm.nih.gov/articles/instance/7643870/bin/NIHMS1581077-supplement-Supp_Combined.pdf |
| SPLiT-seq **Table S12 (all oligos)** | https://pmc.ncbi.nlm.nih.gov/articles/instance/7643870/bin/NIHMS1581077-supplement-Table_S12.xlsx |
| SPLiT-seq Tables S4–S7, S9, Table 1 (saved files S4–S7 are captcha HTML, not xlsx) | https://pmc.ncbi.nlm.nih.gov/articles/instance/7643870/bin/NIHMS1581077-supplement-TableS4.xlsx (etc.) |
| published microSPLiT supplement and **Table S3 (oligos)** | https://pmc.ncbi.nlm.nih.gov/articles/instance/8269303/bin/NIHMS1717416-supplement-Table_S3.xlsx and `.../nihms-1717416.pdf` supplement |

So: **microSPLiT (preprint version) is documented from primary sources; published
SPLiT-seq rests on the upstream page**, cross-checked wherever its oligos are shared with
the microSPLiT table (most of them are).

---

## 1. What it is

No compartments at all: **fixed, permeabilised cells (or nuclei) are the compartments**.
Cells are split into a plate, a well-specific barcode is put on the cDNA inside them, the
cells are pooled, split again, and so on. A transcript's cell of origin is the
*combination* of well barcodes it collected. 🟢 (SPLiT-seq main text)

Four barcodes 🟢: round 1 = **barcoded RT primer** (in-cell RT; the round-1 well also serves
as the sample index), round 2 and round 3 = **in-cell ligation** of barcoded adaptors (the
round-3 oligo also carries the **UMI**), round 4 = the **sublibrary PCR
index** (i7). 96 × 96 × 96 × 24 = 21,233,664 combinations in the paper's arithmetic;
the brain experiment used 48 × 96 × 96 × 14. 🟢 (main text). The round-3 oligo's **5' biotin**
is not in the main text; it comes from the oligo tables (`/5Biosg/`: 🟢 microSPLiT Table S3,
🟡 upstream for SPLiT-seq).

| | New thing here | Where else it turns up |
|---|---|---|
| 1 | **Split-pool ligation barcoding in fixed cells** — no droplets, wells-per-cell or instruments | sci-RNA-seq (combinatorial indexing by RT + PCR), Parse Biosciences (commercial SPLiT-seq) |
| 2 | **Splint (linker) ligation**: a short bridging oligo holds barcode oligo and cDNA 5' end together for T4 DNA ligase, then a **blocking strand** neutralises unused linkers before pooling | the same bridge trick as in [small-RNA ligation](../ref/concepts/small-rna-ligation.md) splints, but DNA–DNA |
| 3 | **Mixed oligo-dT + random-hexamer RT** in every round-1 well | see [reverse transcription](../ref/concepts/reverse-transcription.md) |
| 4 | **Template switching after lysis, on streptavidin beads**, not during the RT | ordinary [template switching](../ref/concepts/template-switching.md), SMART handle |
| 5 | microSPLiT: **in-cell poly(A) polymerase** tailing of bacterial mRNA, lysozyme permeabilisation | PETRI-seq (concurrent, bacterial) |

## 2. Oligos

### 2a. microSPLiT (preprint) — 🟢 verbatim from Auxiliary Table S3 (`869248_media-2.xlsx`)

IDT notation as given in the table: `/5Phos/` (round 1 written `\5Phos\`) = 5' phosphate,
`/5Biosg/` = 5' biotin, `rG` = ribo-G, `+G` = LNA-G. Plates: one example well each; all 96
wells follow the pattern (checked by parsing, 🟡 computed).

```
Round1_01 (A1, dt(15)VN)        \5Phos\ACTGTGGACTCGTAATTTTTTTTTTTTTTTVN
Round1_49 (E1, random hexamer)  \5Phos\ACTGTGGCTGCTTTGTNNNNNN
Round2_01 (A1)                  /5Phos/CATCGGCGTACGACTAACGTGATATCCACGTGCTTGAG
Round3_01 (A1)                  /5Biosg/CAGACGTGTGCTCTTCCGATCTNNNNNNNNNNAACGTGATGTGGCCGATGTTTCG

BC_0335  Round 2 barcode linker          CCACAGTCTCAAGCACG
BC_0340  Round 2 blocking strand         CGTGCTTGAGACTGTGG
BC_0284  Round 3 barcode linker          TACGCCGATGCGAAACATCG
BC_0066  Round 3 blocking strand         GTGGCCGATGTTTCGCATCGGCGTACGACT
BC_0127  Template switching primer       AAGCAGTGGTATCAACGCAGAGTGAATrGrG+G
BC_0243  Adaptor duplex top strand       ACACTCTTTCCCTACACGACGCTCTTCCGATCT
BC_0244  Adaptor duplex bottom strand    GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
BC_0027  PCR primer for sublibrary amplification
         AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
BC_0076  PCR primer (TSBC07), sublibrary index #1 (used with BC_0027)
         CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
  ...    BC_0077-BC_0083 (TSBC08-TSBC14): same, with i7 TCAAGT CTGATC AAGCTA GTAGCC TACAAG TTGACT GGAACT
```

The table's note: random-hexamer and dT(15)VN primers are combined per RT well, "well E1
mixed with A1, etc." 🟢 — so 48 RT wells, each with one dT and one hexamer primer.

### 2b. SPLiT-seq (published) — 🟡 (upstream), Table S12 not fetched

As upstream writes them; it names Table S12 as its source. `[...]` placeholders are upstream's.

```
Oligo-dTVN (Round1_01-48)   /5Phos/ AGGCCAGAGCATTCG[8-bp Round1 barcode]TTTTTTTTTTTTTTTVN
Oligo-randN (Round1_49-96)  /5Phos/ AGGCCAGAGCATTCG[8-bp Round1 barcode]NNNNNN
Round2 (Round2_01-96)       /5Phos/ CATCGGCGTACGACT[8-bp Round2 barcode]ATCCACGTGCTTGAG
Round3 (Round3_01-96)       /5Biosg/ CAGACGTGTGCTCTTCCGATCT[10-bp UMI][8-bp Round3 barcode]GTGGCCGATGTTTCG
BC_0215  Round2 barcode linker      CGAATGCTCTGGCCTCTCAAGCACGTGGAT
BC_0060  Round3 barcode linker      AGTCGTACGCCGATGCGAAACATCGGCCAC
BC_0216  Round2 blocking strand     ATCCACGTGCTTGAGCGCGCTGCATACTTG
BC_0066  Round3 blocking strand     GTGGCCGATGTTTCGCATCGGCGTACGACT
BC_0127  TSO                        AAGCAGTGGTATCAACGCAGAGTGAATrGrG+G
BC_0062  cDNA amplification primer 1   CAGACGTGTGCTCTTCCGATCT
BC_0108  cDNA amplification primer 2   AAGCAGTGGTATCAACGCAGAGT
BC_0076-BC_0083  Indexed library PCR primer 1   CAAGCAGAAGACGGCATACGAGAT[6-bp i7]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
BC_0118  Library PCR primer 2
         AATGATACGGCGACCACCGAGATCTACACTAGATCGCTCGTCGGCAGCGTCAGATGTGTATAAGAGACAG
```

**Checked against the microSPLiT table** 🟡 (computed): the round-2 and round-3 barcode
oligo frames, BC_0066, BC_0127 and the BC_0076–0083 frame are **identical** in both, and
BC_0060 is exactly the reverse complement of the microSPLiT table's BC_0066. What only
upstream attests: the round-1 handle `AGGCCAGAGCATTCG`, BC_0215, BC_0216, BC_0062,
BC_0108, BC_0118. (BC_0062's sequence is, however, the 🟢 5' handle of every round-3
oligo, and BC_0108 the first 23 nt of 🟢 BC_0127.)

### 2c. How they interlock — 🟡 (computed with `lib/`)

**Round-2 adaptor = barcode oligo + linker (splint).** Each linker is two arms:

- SPLiT-seq BC_0215 (30 nt) = `revcomp(AGGCCAGAGCATTCG)` (pairs with the round-1 handle at
  the cDNA's 5' end) + `revcomp(ATCCACGTGCTTGAG)` (pairs with the round-2 oligo's 3' end).
- microSPLiT BC_0335 (17 nt) = `revcomp(ACTGTGG)` (7 nt, the shorter microSPLiT round-1
  handle) + `revcomp` of the round-2 oligo's last **10** nt.
- Either way the round-2 oligo's 3'-OH abuts the cDNA's 5'-phosphate (from the
  `/5Phos/` RT primer): a nick that T4 DNA ligase seals. Result, first strand 5'→3':
  `R2 oligo · R1 primer · cDNA(antisense)`.

**Round-3 adaptor**, same logic:

- SPLiT-seq BC_0060 (30 nt) = `revcomp(CATCGGCGTACGACT)` (round-2 oligo's 5'-phosphate
  end) + `revcomp(GTGGCCGATGTTTCG)` (round-3 oligo's 3' end).
- microSPLiT BC_0284 (20 nt) = the middle of BC_0060 (`BC_0060[5:25]`): 10-nt arms each side.

**Blocking strands.**

- BC_0066 = **exact reverse complement of BC_0060** (30 nt), so it covers the microSPLiT
  BC_0284 completely too.
- microSPLiT BC_0340 = **exact reverse complement of BC_0335** (17 nt).
- SPLiT-seq BC_0216 is **not** the reverse complement of BC_0215: its first 15 nt equal
  the round-2 oligo's 3' handle `ATCCACGTGCTTGAG` (so it pairs with the linker's
  round-2-oligo arm), and its last 15 nt `CGCGCTGCATACTTG` match nothing else in the
  protocol on either strand (checked against the Illumina/Nextera sequences in `lib/` too). 🔴 why.

**Handles that become primer sites.**

- Round-3 5' handle `CAGACGTGTGCTCTTCCGATCT` = BC_0062 = the last 22 nt of
  `illumina.TRUSEQ_READ2`. So the round-3 oligo **is** the Read-2 side of the library; the
  index primers BC_0076–0083 = `illumina.P7` + 6-nt i7 + `illumina.TRUSEQ_READ2` (64 nt)
  extend it.
- BC_0127 (TSO) = `rt.SMART_HANDLE` + `GAAT` + `rGrG+G` (the SMART-seq2 LNA tail,
  `rt.TSO_G_TAIL_LNA`); BC_0108 = `rt.SMART_HANDLE`. Same handle as SMART-seq2's ISPCR.
- SPLiT-seq BC_0118 (70 nt) = `illumina.P5` + `TAGATCGC` + `nextera.READ1_PRIMER`
  (s5 + ME): a single, **fixed i5** — the sublibraries are told apart only by i7.
- microSPLiT BC_0243 = `illumina.TRUSEQ_READ1`; BC_0244 = `revcomp(BC_0243)` minus its
  first base (= `illumina.INDEX2_PRIMER_RC[1:]`), so the duplex has a single **3'-T
  overhang** on BC_0243 — a T-tailed, fully double-stranded (not Y) adaptor. BC_0244 has
  no 5' phosphate in the table. BC_0027 (58 nt) = `illumina.P5` + `TRUSEQ_READ1[4:]`
  (the two share `ACAC`) = `illumina.NEBNEXT_UNIVERSAL_PRIMER`.

**Index orientation.** The i7 is written in the primer as the reverse complement of what
the Index-1 read reports: TSBC07…TSBC14 read as `CAGATC ACTTGA GATCAG TAGCTT GGCTAC
CTTGTA AGTCAA AGTTCC`. The `TSBC07–14` names and these reads look like TruSeq LT indices
7–14; not checked against a list in `lib/`. 🟡

**Round-1 primers (microSPLiT table).** 48 dT primers = `ACTGTGG` + 8-nt barcode +
T15VN (32 nt); 48 hexamer primers = `ACTGTGG` + 8-nt barcode + **`T`** + N6 (22 nt). The
two sets share **no** 8-mer, and neither shares one with the round-2/3 set; the round-2 and
round-3 plates use the **same 96 8-mers** in the same well order. 🟡 (computed from the table)

## 3. Step by step

microSPLiT conditions 🟢 from the preprint Supplementary Materials (and identical where
the published "Methods Summary" repeats them). SPLiT-seq conditions are in its gated
supplement 🔴; microSPLiT defers to them for ligation-plate preparation, cDNA purification
and PCR ("performed according to the SPLiT-seq protocol").

**microSPLiT front end (bacteria)**

1. **Fix**: 4 % formaldehyde in PBS, 4 °C overnight; wash in 100 mM Tris-HCl + RNase inhibitor.
2. **Permeabilise**: 0.04 % Tween-20 in PBS, 3 min on ice; **lysozyme** 2.5 mg/mL (Tris
   pH 7, 50 mM EDTA), 37 °C exactly 15 min.
3. **In-cell polyadenylation**: *E. coli* poly(A) polymerase I (NEB), 10 µL of 10 mM ATP
   in a 100 µL reaction (1 mM final, 🟡 computed), 37 °C 30 min. This is what lets oligo-dT see bacterial mRNA.
4. Vortex, double-filter (10 µm then 1 µm pluriStrainer) — against aggregates.

**Barcoding (both)**

5. **Round 1, in-cell RT**, up to 48 wells, each with its dT + hexamer barcoded primers
   (microSPLiT: 2.5 µM hexamer + 5 µM dT, 1:2), Maxima H Minus RT 20 U/µL, 7.5 % PEG8000,
   500 µM dNTPs, RNase inhibitors; 20 µL, 23 °C 10 min then 50 °C 50 min. 🟢 (microSPLiT)
   First strand 5'→3': `5'P · R1 handle · bc1 · T15VN (or T+N6) · cDNA(antisense)`. The
   round-1 well = sample. Pool; add Triton X-100; spin; vortex, filter, brief sonication
   (10 % power, 5 s) 🟢.
6. **Round 2, in-cell ligation**: split into 96 wells pre-loaded with annealed round-2
   barcode oligo + linker; ligation mix of T4 DNA ligase (NEB), 10X T4 ligase buffer (NEB),
   50 % PEG8000 and RNase inhibitors, combined with the cells 🟢 (final 1× in the wells is
   🟡 inferred) (microSPLiT mix; well volumes and times 🔴, "as in
   SPLiT-seq"). Then the **round-2 blocking strand** is added, cells pooled. 🟡 (upstream;
   the order is in the SPLiT-seq preprint text for its own design 🟢)
7. **Round 3, in-cell ligation** of the biotinylated round-3 oligo (UMI + bc3) with its
   linker, then BC_0066 to block. 🟡 (upstream)
   First strand now: `biotin · R3 handle · UMI · bc3 · R3 3' arm · R2 oligo · R1 primer · cDNA`.
8. **Sublibraries**: wash, count, aliquot; freeze −80 °C overnight (microSPLiT); lyse in
   SDS/EDTA/NaCl + proteinase K, **55 °C 2 h**, which also reverses the crosslinks 🟢
   (microSPLiT; upstream says the same 55 °C for SPLiT-seq).
9. **Capture on streptavidin beads** by the round-3 biotin 🟢 (microSPLiT, "according to
   the SPLiT-seq protocol"). Only cDNA that received the round-3 adaptor is captured — and
   any free round-3 oligo with it. 🟡
10. **Template switch on beads**: Maxima H Minus RT + BC_0127, 50 % PEG8000 (33 µL in
    ~220 µL, ≈7.5 % final 🟡 computed), dNTPs; room temperature 30 min, then 42 °C 90 min
    🟢 (microSPLiT). Whether the non-templated `CCC`
    is added during the in-cell RT (step 5) or during this second RT is not stated. 🔴
    After it, the cDNA's 3' end carries `revcomp(SMART handle + GAAT + GGG)`. 🟡 —
    [template switching](../ref/concepts/template-switching.md)
11. **cDNA PCR on beads**, BC_0062 + BC_0108 (round-3 handle and SMART handle) 🟡
    (upstream; microSPLiT says only "according to the SPLiT-seq protocol" and its table
    lists neither primer 🔴), followed by qPCR-monitored amplification.

**Library: SPLiT-seq** 🟡 (upstream)

12. **Nextera XT tagmentation** of the amplified cDNA. Of the five fragment types upstream
    enumerates, only `s5-ME · insert · …round-3 handle` amplifies with the primers used
    (BC_0118 recognises s5, the index primer recognises the round-3 handle). See
    [Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md).
13. **Sublibrary PCR**: BC_0118 + one of BC_0076–0083 (the i7 = barcode 4).

**Library: microSPLiT preprint** 🟢 (Supplementary Materials, "Fragmentation")

- **Enzymatic fragmentation**: 110 ng cDNA, Enzymatics 5× WGS Fragmentation Mix, 32 °C
  10 min then 65 °C 30 min; double-sided SPRI 0.6–0.8×.
- **Adaptor ligation**: pre-annealed BC_0243/BC_0244 duplex, Enzymatics WGS ligase,
  rapid ligation buffer, 20 °C 15 min; SPRI 0.8×. The 3'-T overhang implies the WGS mix
  end-repairs and dA-tails the fragments 🟡 (kit chemistry; not stated). Every fragment
  end gets the same TruSeq-Read-1 adaptor.
- **Sublibrary PCR**: KAPA HiFi + EvaGreen, BC_0027 + one of BC_0076–0083; 95 °C 3 min;
  cycles of 98 °C 20 s, 67 °C 20 s, 72 °C 3 min, stopped near qPCR saturation; 72 °C 5 min;
  SPRI 0.5–0.7×. Only `adaptor · insert · …round-3 handle` fragments get both primer
  sites. 🟡

## 4. Final libraries — 🟡 (assembled from the oligos above; top strand, 5'→3')

**SPLiT-seq** (identical, base for base in every fixed segment, to upstream's drawn final library — checked by string comparison):

```
5'- P5 · TAGATCGC (fixed i5) · s5 · ME · <cDNA insert, sense> · (A)n · <bc1>
    · CGAATGCTCTGGCCTCTCAAGCACGTGGAT (= BC_0215) · <bc2>
    · AGTCGTACGCCGATGCGAAACATCGGCCAC (= BC_0060) · <bc3> · <UMI 10>
    · AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC (= revcomp TruSeq Read 2) · <i7'> · P7' -3'
```

Fixed bases: 70 before the insert, 158 from bc1 to the end. 🟡 (computed)

**microSPLiT, preprint chemistry** (upstream's microSPLiT library equals this from bc1 onward, but has the Nextera 5' end of SPLiT-seq instead of TruSeq Read 1 — see §6):

```
5'- P5 · TruSeq Read 1 (BC_0027 / BC_0243) · <cDNA insert, sense> · (A)n · <bc1>
    · CCACAGTCTCAAGCACGTGGAT (= revcomp(ACTGTGG) + revcomp(R2 3' arm), 22 nt) · <bc2>
    · AGTCGTACGCCGATGCGAAACATCGGCCAC · <bc3> · <UMI 10>
    · AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC · <i7'> · P7' -3'
```

Fixed bases: 58 before the insert, 150 from bc1 to the end. 🟡 (computed). The barcodes on
the top strand are the reverse complements of the plate sequences; the insert is in mRNA
sense for both dT- and hexamer-primed molecules.

## 5. Sequencing

- **SPLiT-seq** 🟡 (upstream): Read 1 = `nextera.READ1_PRIMER`, 66 cycles, cDNA; Index 1
  = `illumina.INDEX1_PRIMER`, 6 cycles, i7; Read 2 = `illumina.TRUSEQ_READ2`, **94 cycles**.
  No i5 read needed (fixed i5).
- **microSPLiT** 🟢: MiSeq or NextSeq, 150-cycle kits, paired end; **Read 1 74 nt** =
  transcript, **Read 2 86 nt** = UMI + barcodes, **index 6 nt** = sublibrary (barcode 4).
  Read 1 would be the standard TruSeq Read 1 primer given BC_0243 🟡.
- **Read 2 layout** 🟡 (computed): Read 2 is primed on the round-3 handle and reads the
  barcodes **as written in the plate tables**: Read 2 = round-3 oligo minus its 22-nt
  handle, then the round-2 oligo, then the round-1 primer. Positions: UMI 1–10, bc3 11–18,
  bc2 49–56, bc1 **87–94** (SPLiT-seq) or **79–86** (microSPLiT, shorter round-1 handle).
  So 94 and 86 cycles are exactly enough to reach the end of bc1 — both match the stated
  read lengths.

## 6. Upstream vs primary sources

Agreements 🟡 (computed): every microSPLiT oligo upstream lists that also appears in the
preprint Table S3 (round-1/2/3 frames, BC_0335, BC_0284, BC_0340, BC_0066, BC_0127,
BC_0076–0083 frame) is identical; upstream's statement that only the round-2 linker region
is shorter in microSPLiT holds (22 vs 30 nt in the library; the round-3 region is 30 nt in
both); the 55 °C lysis/reverse crosslink matches.

Disagreements and gaps:

- **microSPLiT library prep.** Upstream draws microSPLiT with **Nextera XT** and BC_0118
  (P5-i5-Nextera), BC_0062 and BC_0108, "taken from Supplementary Table S3". The preprint's
  Table S3 lists **none** of those; it lists BC_0243/0244 + BC_0027 (TruSeq adaptor
  ligation after Enzymatics fragmentation), and the preprint methods describe exactly that.
  The published Table S3 and supplement were not fetched, so whether the published
  version switched to Nextera is 🔴.
- **microSPLiT random-hexamer primer**: the table has a constant `T` between barcode and
  N6 (22 nt); upstream writes `ACTGTGG[8-bp]NNNNNN`. 🟢 table.
- Upstream's step 9 names "BC_0018" for the P5 primer it lists as BC_0118 above — a typo
  in one of the two places. 🟡
- SPLiT-seq RT/ligation/PCR conditions and Table S12 itself 🔴 (gated).

## 7. Versions

- **SPLiT-seq preprint (2017)**: RT unbarcoded, then **three** split-ligate-pool rounds
  (96³ = 884,736 combinations), blocking strand after each ligation 🟢 (`105163_v1`). Upstream
  archives a separate page for it; its oligos were not compared here.
- **SPLiT-seq published (2018)**: round 1 moved into barcoded RT primers (dT + hexamer),
  two ligation rounds, i7 sublibrary index as round 4.
- **microSPLiT (2019 preprint / 2021)**: SPLiT-seq barcoding with a 7-nt round-1 handle and
  short linkers, plus bacterial fixation/lysozyme/PAP front end and anti-aggregation
  steps; library by fragmentation + adaptor ligation in the preprint.

## 8. Open questions

- 🔴 SPLiT-seq Table S12 and methods: verify the round-1 handle, BC_0215, BC_0216, BC_0118,
  BC_0062/0108 and all reaction conditions.
- 🔴 Why BC_0216's 3' half (`CGCGCTGCATACTTG`) is unrelated to the round-1 handle while
  every other blocker is the full linker complement — or whether upstream mis-copied it.
- 🔴 Published microSPLiT library prep (Nextera, as upstream, or adaptor ligation, as the
  preprint).
- 🔴 When the template-switch `CCC` is added (in-cell RT vs on-bead RT).
- 🟡 "24 PCR reactions" in the arithmetic vs 8 index primers listed in the microSPLiT table.
- 🟡 Identity of TSBC07–14 with TruSeq LT indices 7–14 — unchecked.

## 9. How this note was made (tool evaluation)

`tools/get_sources.py` saved the upstream page, the SPLiT-seq PMC text and the microSPLiT
preprint with supplements; it reported four SPLiT-seq supplement files as **saved** that are
in fact reCAPTCHA HTML pages (Tables S4–S7), and the microSPLiT PMC page was also a captcha
page. The published microSPLiT text was obtained instead from NCBI efetch (XML), whose
`.txt` twin had to be made by hand because `doctext.py` copied XML through unchanged.
The SPLiT-seq preprint PDF was fetched by hand. `scrape_primers.py` found the microSPLiT
Table S3 first; its `--find` did not locate oligos containing a run of `N` (the round-3
UMI), which were confirmed by grep instead.
