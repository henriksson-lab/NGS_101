# scNOMe-seq / scCOOL-seq — GpC methylase footprinting plus bisulfite sequencing in one cell

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published. Relationships marked 🟡 *(computed)* were worked out with `lib/`; the
> schematic's construct identities are encoded in its protocol module and checked by its
> self-test.

Two independent papers, one principle, two different library chemistries:

**scNOMe-seq** — Pott S. "Simultaneous measurement of chromatin accessibility, DNA
methylation, and nucleosome phasing in single cells." *eLife* 6:e23203 (2017).
doi:[10.7554/eLife.23203](https://doi.org/10.7554/eLife.23203) · PMC5487215 · PMID 28653622.
Single author; human GM12878 and K562 nuclei.

**scCOOL-seq** — Guo F, Li L, Li J, *et al.*, Tang F. "Single-cell multi-omics
sequencing of mouse early embryos and embryonic stem cells." *Cell Research*
27:967–988 (2017). doi:[10.1038/cr.2017.82](https://doi.org/10.1038/cr.2017.82) ·
PMC5539349 · PMID 28621329. Mouse ES cells, oocytes, preimplantation blastomeres.

Built on: **NOMe-seq** (Kelly TK *et al.* 2012, bulk; GpC methyltransferase M.CviPI +
bisulfite sequencing) and, for scCOOL-seq, **scBS-seq / PBAT**
(Smallwood *et al.* 2014, see [scBS-seq](../scbs-seq__10.1038+nmeth.3035/01_scbs-seq.md)),
which scCOOL-seq cites as its ref. 16. The same M.CviPI step is later combined with
scM&T-seq into scNMT-seq (separate catalogue entry `scnmt-seq__10.1038+s41467-018-03149-4`).

Sources read (in `_data/sources/scnome-seq-sccool-seq__10.7554+eLife.23203/`, never committed):

| File | What | Used for |
|---|---|---|
| `eLife.23203_PMC5487215.html.txt` | scNOMe-seq full text (PMC) | **all scNOMe-seq methods** ("Materials and methods", line 1523 ff.) |
| `elife-23203-supp1-v1.xlsx.txt` | eLife Supplementary file 1 (fetched by hand from `cdn.elifesciences.org`) | **every scNOMe-seq oligo** (sheet "Supplemental Table 2"); per-cell barcode, read type, sequencer (sheet "Supplemental Table 1") |
| `cr201782.html.txt` | scCOOL-seq full text (nature.com HTML, fetched by hand) | **all scCOOL-seq methods and its two random primers** ("Single-cell COOL-seq library preparation and sequencing", line 934) |
| `upstream_scNOMe_scCOOL.html.txt` | scg_lib_structs page | second account: principle only; defers to its scBS-seq page for the library |
| `cr.2017.82_PMC5539349.html`, `.xml` | PMC copy of the scCOOL-seq paper (PMC5539349; arrived late from `get_sources.py`, no `.txt` twin) | not used — duplicate of `cr201782.html` |
| `eLife.23203_PMC5487215.xml`, `eLife.23203_PMC5487215_supplementary.zip` | JATS duplicate of the eLife text; PMC supplement bundle | not used (supplement read from the `.xlsx`) |

Not fetched / not available:

| What | Why it matters |
|---|---|
| Zymo **Pico Methyl-Seq Library Prep Kit** manual / adapter sequences | the scNOMe-seq random primers, PrepAmp and amplification primers are all **proprietary kit reagents** — no sequence in either paper 🔴 |
| scCOOL-seq supplementary PDFs (`…MOESM357–361_ESM.pdf`, linked from https://www.nature.com/articles/cr201782) | figures and motif tables; the methods text gives no further oligos, so probably not chemistry 🟡 |
| NEB "Illumina Reverse indexed primer" used by scCOOL-seq | the paper names only the vendor; the sequence is not given 🔴 |
| Kelly *et al.* 2012 bulk NOMe-seq (doi:10.1101/gr.143008.112) | parent method; not needed for the library chemistry |

---

## 1. What it is

Mammalian DNA carries endogenous methylation at **CpG** but essentially none at **GpC**.
Both methods treat permeabilised nuclei with the bacterial GpC methyltransferase
**M.CviPI** + SAM 🟢: only GpCs that are **not wrapped in a nucleosome or covered by a
bound factor** get methylated. Bisulfite then converts every unmethylated C to U, and a
single whole-genome bisulfite library reads **two layers at once** 🟢:

- **WCG / HCG** methylation = endogenous DNA methylation (scCOOL-seq says WCG, scNOMe-seq
  CpG minus GCG);
- **GCH** methylation = accessibility (M.CviPI footprint), along the whole read, which also
  gives nucleosome positioning / phasing and TF footprints (CTCF in scNOMe-seq).
- **GCG** is ambiguous (both GpC and CpG) and is dropped by both papers 🟢; scCOOL-seq also
  drops CCG 🟢 (off-target M.CviPI activity).

scCOOL-seq adds two read-count analyses from the same data — CNV (HMMcopy) and **ploidy
from a λ-DNA spike-in** (≈1 pg per cell) 🟢 — which is what the upstream page means by
"extra analysis to investigate CNV and ploidy".

| | New thing here | Consequence |
|---|---|---|
| 1 | **Enzymatic accessibility mark written into the DNA itself** (5mC at GpC), read by bisulfite — no Tn5, no DNase cut | accessibility and CpG methylation come from the **same read**; coverage does not depend on accessibility, so "closed" and "missing" can be told apart 🟢 (scNOMe-seq intro) |
| 2 | scNOMe-seq: M.CviPI on **bulk nuclei, then FACS** one nucleus per well; library from a **commercial kit** (Zymo Pico Methyl-Seq) | no custom chemistry except the index PCR primers |
| 3 | scCOOL-seq: lysis, M.CviPI, protease and bisulfite **in one tube per cell**; library by **scBS-seq-style PBAT** with a **different second-strand primer** (TruSeq Read 2 handle, not scBS-seq oligo2) | an off-the-shelf NEB indexed primer can be used for PCR 🟡 (inferred from the handle sequence) |

There is **no cell barcode and no UMI** in either chemistry: one cell per well / tube,
identified by the sample index of the final PCR primer. None of the repository concept
pages (Tn5, template switching, ligation, RT, padlock) applies; the library chemistry is
DNA-polymerase random priming on bisulfite-converted single-stranded DNA, as in scBS-seq. 🟡

## 2. Oligos

### scNOMe-seq — 🟢 verbatim from Supplementary file 1, sheet "Supplemental Table 2"

No modifications are written. The paper says "primers were ordered from IDT".

```
forward_p5   AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT          58 nt
INDEX16      CAAGCAGAAGACGGCATACGAGATGGACGGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT    64 nt
INDEX19      CAAGCAGAAGACGGCATACGAGATTTTCACGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
INDEX21      CAAGCAGAAGACGGCATACGAGATCGAAACGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
INDEX22      CAAGCAGAAGACGGCATACGAGATCGTACGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
INDEX23      CAAGCAGAAGACGGCATACGAGATCCACTCGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
INDEX25      CAAGCAGAAGACGGCATACGAGATATCAGTGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
INDEX26      CAAGCAGAAGACGGCATACGAGATGCTCATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
INDEX27      CAAGCAGAAGACGGCATACGAGATAGGAATGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
INDEX28      CAAGCAGAAGACGGCATACGAGATCTTTTGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
INDEX29      CAAGCAGAAGACGGCATACGAGATTAGTTGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
INDEX30      CAAGCAGAAGACGGCATACGAGATCCGGTGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
INDEX31      CAAGCAGAAGACGGCATACGAGATATCGTGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
INDEX32      CAAGCAGAAGACGGCATACGAGATTGAGTGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
INDEX33      CAAGCAGAAGACGGCATACGAGATCGCCTGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
INDEX34      CAAGCAGAAGACGGCATACGAGATGCCATGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
INDEX35      CAAGCAGAAGACGGCATACGAGATAAAATGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
INDEX36      CAAGCAGAAGACGGCATACGAGATTGTTGGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
INDEX37      CAAGCAGAAGACGGCATACGAGATATTCCGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
INDEX38      CAAGCAGAAGACGGCATACGAGATAGCTAGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
INDEX39      CAAGCAGAAGACGGCATACGAGATGTATAGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
INDEX40      CAAGCAGAAGACGGCATACGAGATTCTGAGGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
INDEX43      CAAGCAGAAGACGGCATACGAGATGCTGTAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
INDEX44      CAAGCAGAAGACGGCATACGAGATATTATAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
INDEX47      CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
```

Kit-internal oligos (Pico Methyl-Seq random primers, PrepAmp and "first amplification"
primers): **not published** 🔴.

### scCOOL-seq — 🟢 verbatim from the Cell Research methods (`cr201782.html.txt` line 934)

```
random primer 1 (first strand)   5′-Biotin-CTACACGACGCTCTTCCGATCTNNNNNNNNN-3′     31 nt
random primer 2 (second strand)  5′-AGACGTGTGCTCTTCCGATCTNNNNNNNNN-3′         30 nt
PCR forward                      "Illumina Forward PE1.0 primer" — sequence not given 🔴
PCR reverse                      "Illumina Reverse indexed primer (New England Biolabs)" — sequence not given 🔴
```

The paper writes the biotin as `5′-Biotin-` 🟢; scBS-seq writes `[Btn]`, and the upstream
scBS-seq page `/Bio/`. The "PE1.0" forward primer is presumably the 58-nt scBS-seq PE1.0
(= scNOMe-seq forward_p5, below) 🟡 — the name matches, the sequence is not printed here. Primer 1 bases are identical to scBS-seq oligo1. 🟡 (computed: same 31-mer)

### Upstream vs. papers

The upstream page draws no oligos of its own for this entry: it states the principle and
says scNOMe-seq uses the Pico Methyl-Seq kit while scCOOL-seq "simply uses scBS-seq",
and that the final libraries of the two are the same.

| Point | Papers | Upstream | Verdict |
|---|---|---|---|
| M.CviPI marks accessible GpC; CpG = endogenous methylation | both papers | same | agree 🟢 |
| scNOMe-seq library = Zymo Pico Methyl-Seq kit | eLife methods | same | agree 🟢 |
| scCOOL-seq library = scBS-seq | "single-cell PBAT strategy, as previously described" (cites Smallwood 2014), **but** the second-strand primer is `AGACGTGTGCTCTTCCGATCT`-N9, not scBS-seq oligo2 `TGCTGAACCGCTCTTCCGATCT`-N9, and the reverse PCR primer is an NEB indexed primer, not iPCRTag | "simply uses scBS-seq" | **partial disagreement**: same strategy, different read-2 handle and index primer 🟡 (computed from the two sequences) |
| final libraries of scNOMe-seq and scCOOL-seq "the same" | scNOMe-seq index primers end in the TruSeq Read 2 arm; scCOOL-seq primer 2 is the 3′ 21 nt of that arm; both forward primers sit on the TruSeq Read 1 site | same | **agree** for the read-2 side as far as can be checked; but the scNOMe-seq insert-proximal ends come from kit primers 🔴 |
| CNV / ploidy analysis in scCOOL-seq | λ-DNA spike-in, HMMcopy | "extra analysis" | agree; upstream does not mention the λ spike-in |
| number of oligo-1 priming rounds in scCOOL-seq | not stated (one annealing/extension is described) | inherits scBS-seq's five | **unresolved** 🔴 |

### How the oligos interlock — 🟡 (computed with `lib/illumina`)

- **forward_p5 = `illumina.P5` + `illumina.TRUSEQ_READ1[4:]`** — identical to
  `illumina.NEBNEXT_UNIVERSAL_PRIMER` and to scBS-seq's PE1.0. Its 3′ 22 nt,
  `CTACACGACGCTCTTCCGATCT`, are exactly the handle of scCOOL-seq random primer 1 (= the
  3′ 22 nt of `TRUSEQ_READ1`). So scCOOL-seq's "Illumina Forward PE1.0 primer" reaches
  primer 1 with 22 nt of overlap, as in scBS-seq.
- **Every INDEXnn = `illumina.P7` (24 nt) + 6-nt index + `illumina.TRUSEQ_READ2`
  (34 nt, `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`)** = 64 nt, for all 24 primers. The
  reverse complement of that 34-nt arm, minus its first base, is `illumina.INDEX1_PRIMER`
  — so the index read is primed off the standard TruSeq index-1 primer.
- **Index orientation.** The 6 nt written in the primer is the **reverse complement of the
  i7 read**: e.g. INDEX16 is written `GGACGG`, so the index read is `CCGTCC`. These
  i7 reads (`CCGTCC`, `GTGAAA`, `GTTTCG`, `CGTACG`, …) look like Illumina's TruSeq LT
  index set with the same numbers; `lib/` has no TruSeq 6-nt index table, so this was not
  checked mechanically. 🟡 The 24 written indices have a minimum pairwise Hamming
  distance of 3. 🟡 (computed)
- **scCOOL-seq primer 2 handle** `AGACGTGTGCTCTTCCGATCT` (21 nt) is
  `illumina.TRUSEQ_READ2[13:]` — the 3′ 21 nt of the same TruSeq Read 2 arm that ends
  every scNOMe-seq INDEX primer. A standard NEB/TruSeq-style indexed P7 primer
  (P7 + i7 + `TRUSEQ_READ2`) therefore anneals to it with 21 nt of overlap; that is
  presumably why the authors could buy the reverse primer from NEB. 🟡
- Primer 1 and primer 2 handles share their 3′ 13 nt `GCTCTTCCGATCT` (the TruSeq Y-stem
  region; `illumina.STEM_COMPLEMENT` + `T`). scBS-seq's oligo2 shares 14 nt with oligo1
  instead (`CGCTCTTCCGATCT`); the two scCOOL-seq handles differ more. 🟡 (computed)
- Supplementary Table 1 uses 24 distinct barcodes, every one of them present in Table 2;
  within each sequencing pool (a = GM12878, 19 cells; b = K562, 11 cells) no barcode
  repeats, but barcodes 25, 29, 30, 33, 35, 36 are reused across the two pools. 🟡 (computed)

## 3. Step by step

### 3a. scNOMe-seq (eLife) — 🟢 unless marked

1. **Nuclei.** 2–5 × 10⁶ cells, washed in PBS, swollen 10 min on ice in 10 mM Tris-HCl
   pH 7.4, 10 mM NaCl, 3 mM MgCl₂; IGEPAL CA-630 to 0.025 %; 15 Dounce strokes;
   800 × g, 4 °C; two washes in detergent-free lysis buffer.
2. **GpC methylation, in bulk.** 1 × 10⁶ nuclei in 150 µL: 1× GpC MTase buffer (NEB),
   0.32 mM SAM, 50 µL M.CviPI (4 U/µL), 37 °C 8 min; then +25 µL enzyme and +0.7 µL
   32 mM SAM, another 8 min at 37 °C. Stop by diluting with 750 µL PBS and spinning
   (800 × g) rather than by heat inactivation, so the nuclei stay intact for sorting
   (the paper's stated reason is to avoid disrupting the nuclei).
3. **Sort.** Nuclei in PBS + Hoechst 33342; FACS (BD FACSAria) one nucleus per well of a
   96-well plate, gated on scatter for singlets and on Hoechst for **G1 DNA content**.
   Wells pre-filled with 19 µL 1× M-Digestion buffer (Zymo EZ DNA Methylation Direct)
   + 1 mg/mL proteinase K.
4. **Lyse** 50 °C 20 min.
5. **Bisulfite** in the well: +130 µL CT Conversion reagent; 98 °C 8 min, 65 °C 3.5 h.
   Clean-up with the Direct kit, **eluted in only 8 µL**. Unmethylated C → U; both
   M.CviPI-methylated GpC and endogenous mCpG stay C. 🟡 (standard bisulfite chemistry)
6. **Library: Zymo Pico Methyl-Seq kit**, low-input instructions, with two changes:
   random primers **diluted 1:2** before the pre-amplification step, and the **first
   amplification extended to 10 cycles** in total. What the kit's steps add to the
   fragment ends is not published 🔴. That the kit leaves a TruSeq Read 1 tail
   (`…CTACACGACGCTCTTCCGATCT`) and a TruSeq Read 2 tail (`AGATCGGAAGAGCACACGTC…`) on
   the amplified fragments is **inferred** from the fact that forward_p5 and the INDEX
   primers anneal there. 🟡
7. **Barcoded amplification** with forward_p5 + one INDEXnn per cell (whether this is the
   kit's final amplification step with the primers swapped in is not stated 🟡; cycles not
   given 🔴).
   AMPure XP 1:1; Bioanalyzer; only libraries with mean fragment size **> 150 bp**
   pooled. Negative controls: GM12878 nuclei processed **without M.CviPI** (7 of 19
   GM12878 libraries; Supplementary file 1, Supplemental Table 1, "GpC MTase = No").

### 3b. scCOOL-seq (Cell Research) — 🟢 unless marked

1. **Pick** one cell / blastomere by mouth pipette into 3.5 µL ice-cold lysis buffer:
   50 mM Tris-HCl pH 7.4, 50 mM NaCl, 10 mM DTT, 0.25 mM EDTA, 0.25 mM PMSF, 0.5 % NP-40,
   **plus 1 pg λ DNA** (ploidy spike-in). 10 min on ice.
2. **GpC methylation in the same tube**: add M.CviPI + SAM to 5 µL final, 1 U/µL
   M.CviPI, 160 µM SAM; 37 °C 30 min; heat-inactivate 65 °C 20 min. The λ DNA is
   naked and should be methylated at every GpC — an internal positive control for
   M.CviPI 🟡 (inferred; the paper uses λ only for ploidy).
3. **Release DNA**: +0.5 µL 20 mg/mL Qiagen protease, 50 °C 3 h.
4. **Bisulfite**: MethylCode Bisulfite Conversion Kit (Invitrogen), per manufacturer.
5. **First strand**: anneal random primer 1 (5′-biotin, Read 1 handle, N9) and extend
   with **Klenow exo-**. Number of rounds, primer concentration and temperature ramp
   not stated 🔴 (scBS-seq used five rounds). Product:
   `5′-biotin · CTACACGACGCTCTTCCGATCT · N9 · copy of a converted strand -3′`. 🟡
6. **Exonuclease I** digests leftover primer 1; AMPure XP.
7. **Capture** on **M-280 streptavidin** Dynabeads; "the original bisulfite-converted DNA
   templates were removed" (washes not specified 🔴; scBS-seq used 0.1 N NaOH).
8. **Second strand**: random primer 2 (Read 2 handle `AGACGTGTGCTCTTCCGATCT`, N9) +
   Klenow exo- on the beads. Product: duplex with the Read 1 handle at one end and the
   Read 2 handle at the other. 🟡
9. **PCR on the beads, 13 cycles**: Illumina Forward PE1.0 + NEB Illumina reverse indexed
   primer, KAPA HiFi HotStart. Two AMPure XP clean-ups; Fragment Analyzer; qPCR
   quantification; pool.

## 4. Final library — 🟡 (assembled from the oligos above)

**scCOOL-seq** (top strand, 5′→3′; the reverse primer's internal sequence is assumed to be
the standard P7 + i7 + `TRUSEQ_READ2` layout of NEB's indexed primers 🟡):

```
5'- P5 · TruSeq Read 1 site · N9 (primer 1) · <bisulfite-converted insert> · rc N9 (primer 2) · rc TruSeq Read 2 arm · <i7> · P7' -3'
```

By parts:

```
AATGATACGGCGACCACCGAGATCTACAC        P5 (illumina.P5; last 4 nt ACAC = first 4 of Read 1)
TCTTTCCCTACACGACGCTCTTCCGATCT        rest of TRUSEQ_READ1; last 22 nt = primer 1 handle
NNNNNNNNN                            primer 1 random 9-mer (not genomic)
<insert>                             bisulfite-converted DNA, either strand orientation (non-directional)
NNNNNNNNN                            primer 2 random 9-mer, as reverse complement
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC   rc of TRUSEQ_READ2; first 21 nt = rc of primer 2 handle
<i7>                                 sample index (NEB primer; not given)
ATCTCGTATGCCGTCTTCTGCTTG             P7' (illumina.P7_RC)
```

**scNOMe-seq**: the outer parts are fixed by the published primers —

```
5'- P5 · TruSeq Read 1 site · [kit random-primer remnant?] · <bisulfite-converted insert> · [kit remnant?] · rc TruSeq Read 2 arm · <6-nt i7 = rc of written index> · P7' -3'
```

Outer fixed length (P5 + rest of Read 1 + rc Read 2 arm + P7') = 58 + 34 + 24 = 116 nt,
plus the 6-nt index. 🟡 (computed). Whether a random-primer stretch sits between the
Read 1 / Read 2 sites and the insert (as N9 does in PBAT) is not published 🔴; the
paper's 6-nt 5′ clipping "to avoid mismatches introduced by amplification" suggests
priming-derived non-genomic bases at the read start. 🟡

## 5. Sequencing

**scNOMe-seq** 🟢: **2 × 100**, HiSeq 2500 rapid mode (K562, pool b) and HiSeq 4000
(GM12878, pool a); 6-nt i7 index read. Mates mapped **independently as single reads**
(better mapping rate); Trim Galore `--clip_R1 6` (GM12878: 6 bp clipped from either end),
quality 30, min length 20; Bismark **`--non_directional`** against hg38; samtools rmdup;
`bismark_methylation_extractor --ignore 6 --CX`, then coverage2cytosine; GCG positions
removed; reads with > 3 unconverted non-CpG/non-GpC cytosines discarded. 🟢

**scCOOL-seq** 🟢: **2 × 150**, HiSeq 2500. Trim Galore removes the **first 9 bases**
(the N9); Bismark paired-end, **non-directional**, mm9; unmapped pairs re-aligned single-end;
samtools rmdup. WCG = endogenous methylation, GCH = accessibility; GCG and CCG excluded.

Read 1 starts in the primer-1 N9 (scCOOL-seq) and Read 2 in the primer-2 N9, read with the
standard TruSeq Read 1 / Read 2 primers 🟡 — the same layout as scBS-seq apart from the
read-2 side, which here is the TruSeq arm (scBS-seq used a custom read-2 primer).

## 6. Open questions

- **Pico Methyl-Seq kit adapters** — what exactly flanks the insert in scNOMe-seq
  libraries (random-primer length, any tail between Read 1/Read 2 site and insert)? 🔴
- scCOOL-seq: how many rounds of primer-1 extension, at what concentrations, and how the
  template strand was removed after capture (NaOH as in scBS-seq?) — not stated 🔴.
- The NEB "Illumina Reverse indexed primer" set used by scCOOL-seq and its index list 🔴.
- scNOMe-seq index PCR cycle number 🔴.
- The scNOMe-seq i7 set matches Illumina's TruSeq LT numbering by eye of the sequences
  only; a TruSeq 6-nt index table in `lib/illumina` would let this be checked. 🟡
- M.CviPI conditions differ: scNOMe-seq ≈1.3 U/µL (200 U in 150 µL) for 2 × 8 min on
  10⁶ intact nuclei, then a second enzyme addition; scCOOL-seq 1 U/µL for 30 min on one
  NP-40-lysed cell. 🟡 (computed from the stated volumes) Whether the shorter scNOMe-seq
  incubation saturates accessible GpCs is not addressed in the methods read here.
