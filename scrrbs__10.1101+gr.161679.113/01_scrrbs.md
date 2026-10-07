# scRRBS — single-cell reduced-representation bisulfite sequencing

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived, inferred, or only from
> the secondary upstream page · 🔴 not published / not available. Relationships marked
> 🟡 *(computed)* were worked out with `lib/` while writing this note; they are not yet
> asserted in a self-test, because this protocol has no `tools/` module yet.

**scRRBS** — Guo H, Zhu P, Wu X, Li X, Wen L, Tang F. "Single-cell methylome landscapes
of mouse embryonic stem cells and early embryos analyzed using reduced representation
bisulfite sequencing." *Genome Research* 23:2126–2135 (2013).
doi:[10.1101/gr.161679.113](https://doi.org/10.1101/gr.161679.113) · PMID 24179143 ·
PMC3847781. Data: GEO GSE47343.

Step-by-step protocol: Guo H *et al.* (same group; full author list not in the fetched
sources). "Profiling DNA methylome landscapes of mammalian cells with single-cell
reduced-representation bisulfite sequencing." *Nature Protocols* 10:645–659 (2015).
doi:[10.1038/nprot.2015.039](https://doi.org/10.1038/nprot.2015.039) · PMID 25837417.

Sources read (fetched by `tools/get_sources.py`, into
`$CHEM_DATA/sources/scrrbs__10.1101+gr.161679.113/`, never committed):

| File | What | Used for |
|---|---|---|
| `PMC3847781.html.txt` | Genome Research paper, full text (PMC) | **all wet-lab conditions** (section "Construction of single-cell RRBS sequencing libraries") |
| `upstream_scRRBS.html.txt` | scg_lib_structs page (drawn from the *Nature Protocols* paper) | **every oligo sequence**, library structure, read layout — secondary source |
| `nprot.2015.039_…MOESM430_ESM.pdf.txt` | *Nat Protoc* Supplementary Figures 1–5 | nothing chemical (CpG coverage, methylation histograms) |
| `nprot.2015.039_…MOESM429_ESM.zip` → `zip/*.sh`, `SingleC_MetLevel.pl` | *Nat Protoc* analysis scripts | trimming / alignment parameters (§6) |

Not obtained — fetch by hand:

| What | URL | Why it matters |
|---|---|---|
| *Nature Protocols* main text (paywalled) | https://www.nature.com/articles/nprot.2015.039 | the only primary source of the adapter sequences, primer names (QP1/QP2), reagent table, PCR cycle numbers of the 2015 version |
| Genome Research supplement (not linked from PMC) | https://genome.cshlp.org/content/23/12/2126/suppl/DC1 | Supplemental Table 5 (nested-PCR primers for the MspI/HpaII sperm test), Supplemental Figs. 9–11 |

The Genome Research paper itself **gives no oligo sequence**: it says only "Illumina
standard premethylated indexed adaptors" 🟢. Every sequence below therefore rests on the
upstream page alone and is 🟡.

---

## 1. What it is

Bulk RRBS (Meissner/Gu: MspI digest → end repair/A-tail → ligate methylated adapters →
bisulfite → PCR) shrunk to one cell. The change is not new chemistry but **removing every
purification before bisulfite**: lysis, lambda spike-in, MspI digestion, fill-in/A-tailing,
adapter ligation and the bisulfite reaction itself all happen **in the same tube** 🟢; the
first clean-up is the Zymo column after bisulfite, with 10 ng tRNA as carrier 🟢.

No cell barcode, no UMI. One cell = one library, identified by the **6-nt TruSeq index in
the ligated adapter** 🟡 (upstream; the paper only says "indexed"). No reverse
transcription, tagmentation, template switching or circularization is involved, so none
of the `ref/concepts/` pages applies; the closest relative there is
[Tn5 tagmentation](../ref/concepts/tn5-tagmentation.md) only as the contrasting way
other single-cell DNA methods attach adapters.

| | Feature | Note |
|---|---|---|
| 1 | One-tube lysis → digestion → ligation → bisulfite | the paper's stated reason: avoiding DNA loss in purification steps 🟢; a diploid mouse cell holds ~6 pg of DNA 🟡 (general knowledge, not in the source) |
| 2 | **60 fg unmethylated lambda DNA** spiked into every cell | per-cell bisulfite conversion rate: 99.2 % on average (97.7–99.9 %) 🟢 |
| 3 | **Adapters ligated before bisulfite**, so their cytosines must be methylated | otherwise the adapter Cs would become U/T and P5/P7 would no longer match 🟡 |
| 4 | Uracil-tolerant **PfuTurbo Cx** for the first PCR | ordinary proofreading polymerases stall on uracil 🟢 (enzyme) / 🟡 (reason) |

## 2. Oligos

All 🟡 — verbatim from the upstream scg_lib_structs page (which follows *Nature
Protocols*), not checkable against a primary text. Upstream: "cytosines are methylated" in
both adapter strands; `/phos/` = 5' phosphate; `[6-bp index]` = the sample (= cell) index.

```
Universal TruSeq adaptor   AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT
Indexed TruSeq adaptor     /phos/GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[6-bp index]ATCTCGTATGCCGTCTTCTGCTTG

QP1 (Illumina P5)          AATGATACGGCGACCACCGAGATCTACAC
QP2 (Illumina P7)          CAAGCAGAAGACGGCATACGAGAT

TruSeq Read 1 primer       ACACTCTTTCCCTACACGACGCTCTTCCGATCT
TruSeq Read 2 primer       GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
Index read primer          GATCGGAAGAGCACACGTCTGAACTCCAGTCAC
```

The names QP1/QP2 are attributed by upstream to the *Nature Protocols* paper; the
Genome Research paper names no PCR primer at all. 🔴 (primary text not obtained)

### How they interlock — 🟡 (computed)

- **Universal adaptor = `illumina.TRUSEQ_P5_FULL`** exactly (58 nt) = `illumina.P5`
  (29 nt) + `illumina.TRUSEQ_READ1` (33 nt), overlapping by 4 nt (`ACAC`). It ends in the
  3'-**T** overhang that pairs with the A-tail.
- **Indexed adaptor = `illumina.INDEX1_PRIMER` (33 nt) + 6-nt index + `illumina.P7_RC`
  (24 nt)** = 63 nt. Its first 33 nt are the reverse complement of `TRUSEQ_READ2`
  minus its 5'-terminal `A` (revcomp(Read 2 primer) = `A` + indexed adaptor[:33]) — that
  `A` is supplied by the A-tail of the insert.
- **Y-adapter**: the reverse complement of the universal adaptor's 3' end, after its single
  3'-T overhang, matches the indexed adaptor's first **12 nt `GATCGGAAGAGC`**; the rest
  of the two oligos is non-complementary — the standard TruSeq fork.
- The 5' phosphate is on the indexed adaptor only: it ligates to the 3'-A of an insert
  strand; the universal adaptor's 3'-T ligates to the insert's other strand, whose 5'
  phosphate comes from the MspI cut.
- Because the adaptor already contains full P5 and P7, **QP1/QP2 = `illumina.P5` /
  `illumina.P7`** are enough to amplify; no index is added by PCR.
- Read primers are the stock TruSeq ones: Read 1 = `TRUSEQ_READ1`, Read 2 =
  `TRUSEQ_READ2`, index read = `INDEX1_PRIMER`.
- Methylated cytosines to protect: 21 C in the universal adaptor, 17 C in the fixed part
  of the indexed adaptor (plus any C in the index).

## 3. Step by step

Conditions 🟢 from Genome Research Methods unless marked. All steps through bisulfite in
one PCR tube in a thermocycler.

1. **Pick and lyse**: single cell by mouth pipette into 5 µL lysis buffer (20 mM
   Tris-EDTA pH 8.0, 20 mM KCl, 0.3 % Triton X-100) + 1 mg/mL Qiagen protease;
   **60 fg unmethylated lambda DNA** spiked in. 50 °C 3 h, protease killed 75 °C 30 min.
2. **MspI digest**: 9 U MspI (C^CGG), 18 µL, 37 °C 3 h. Every fragment starts `CGG` on
   both strands with a 5'-`CG` overhang. MspI is methylation-insensitive, so cutting does
   not depend on the CpG's state. 🟡
3. **Fill-in + A-tail**: 5 U Klenow exo⁻, 1 mM dATP + 0.1 mM dGTP + 0.1 mM dCTP, 20 µL.
   The 3' recessed end is filled with `CG` and then given an extra `A`; each strand is now
   `5'-CGG…CCGA-3'`. 🟡 (mechanism; drawn the same way upstream) The filled-in C is from
   **unmethylated dCTP**, so it will read as T whatever the genome had there — an
   artificial position at the 3' end of every fragment. 🟡
4. **Ligate** pre-methylated indexed Y-adapters: 30 U high-concentration T4 DNA ligase,
   25 µL. Adapter amount, ligation time and temperature are not given 🔴.
5. **Bisulfite** ("MethyCode" kit as the paper spells it — Invitrogen MethylCode), 150 µL: denature 98 °C 10 min, convert
   64 °C 2.5 h. On-column desulfonation, Zymo-Spin column with **10 ng tRNA carrier**,
   elute in 30 µL. Unmethylated C → U in the insert; the methylated adapter Cs survive. The
   two strands are no longer complementary and separate.
6. **PCR round 1**: 50 µL, 1 U **PfuTurbo Cx** (reads through uracil): 95 °C 2 min;
   25 × (95 °C 20 s, 60 °C 30 s, 72 °C 1 min).
7. **PCR round 2**: 50 µL, 1 U Phusion: 98 °C 2 min; 22 × (98 °C 10 s, 60 °C 30 s,
   72 °C 1 min). Primers for either round are not named in the paper 🔴; upstream draws a
   single P5/P7 (QP1/QP2) PCR 🟡. Whether round 2 is on a dilution or clean-up of round 1
   is not said 🔴.
8. **Size-select** 200–500 bp on a 12 % native polyacrylamide TBE gel; Fragment Analyzer
   QC; qPCR quantification.

In the first PCR cycle only P7 (QP2) can anneal: each converted strand carries the
universal adaptor at its 5' end and the P7 complement at its 3' end; P5 finds a site only
on the copy. Both original strands give the same structure. 🟡 (upstream, consistent with
the adapter design)

Side assay (not part of the library): single sperm lysed the same way, digested with MspI
or the methylation-sensitive **HpaII**, then nested PCR (30 + 30 cycles) — primers in
Supplemental Table 5, not obtained 🔴.

## 4. Final library — 🟡 (assembled from the oligos; agrees with upstream's drawing)

```
5'- P5 · TruSeq Read 1 · <insert: Y-GG … Y-T-G> · A · INDEX1_PRIMER seq (= Read 2') · <i7, 6 nt> · P7' -3'
```

- `P5 · TruSeq Read 1` = the universal adaptor (58 nt); `INDEX1_PRIMER · i7 · P7'` = the
  indexed adaptor (63 nt). Adapters add **121 nt** (+ the A-tail) to the insert.
- Insert: starts at the MspI site, `CGG` if that CpG was methylated, `TGG` if not
  (bisulfite); ends `…C/T` `T` `G`, the `T` being the artificial filled-in C; then the `A`
  tail. Non-CpG C are expected to read as T.
- Upstream's drawing matches this exactly (`…GATCTTGG…TTGAGATCGG…`), showing the
  unmethylated case. Size-selected 200–500 bp total 🟢 → insert roughly 80–380 bp 🟡.

## 5. Sequencing

🟢 Paired-end on an **Illumina HiSeq 2000**; read length not stated in the Genome
Research paper 🔴.

🟡 (upstream, standard TruSeq):

| Read | Primer | Reads |
|---|---|---|
| Read 1 | TruSeq Read 1 | insert from the MspI end: first base = methylation call of that CpG, then `GG` |
| Index (i7) | INDEX1_PRIMER | the 6-nt sample index = **cell identity** |
| Read 2 | TruSeq Read 2 | the other end of the same converted strand (directional library: reads are complementary to the converted strand); position 2 is the filled-in, artificial C |

No i5 read: the universal adaptor has no index.

## 6. Analysis hints from the scripts — 🟢 (MOESM429 zip)

- `trim_galore --quality 20 --phred33 --stringency 3 --length 36 --rrbs --paired --trim1` — the
  `--rrbs` switch exists to remove the filled-in, unmethylated-dCTP positions at the
  fragment ends (§3 step 3) 🟡 (purpose; standard Trim Galore behaviour).
- `bismark` with default parameters on mm9 🟢 (paper); the alignment script calls
  bowtie 1.0.0 (`--path_to_bowtie <bowtie-1.0.0>`) 🟢 (script), although the genome
  preparation script uses `--bowtie2` 🟢 (script) — an inconsistency in the supplied
  scripts. Default Bismark mode is directional 🟡. Lambda genome (48,502 bp) as an extra
  reference for the conversion rate 🟢 (paper).
- Calls: ≥3× coverage per CpG; ≥90 % → methylated, ≤10 % → unmethylated, others
  discarded 🟢 (paper).

## 7. Agreements and disagreements: paper vs upstream

| Point | Genome Research 2013 | Upstream (after *Nat Protoc* 2015) |
|---|---|---|
| Adapters | "Illumina standard premethylated indexed adaptors", no sequence | full TruSeq Y-adapter sequences, Cs methylated, 6-nt index — consistent |
| Enzyme / fill-in / A-tail | MspI, Klenow exo⁻, dATP/dGTP/dCTP | same, drawn |
| PCR | **two rounds** (PfuTurbo Cx 25 cyc, then Phusion 22 cyc), primers unnamed | one P5/P7 (QP1/QP2) PCR drawn, no polymerase, no cycle numbers |
| Lambda spike-in, tRNA carrier, PAGE size selection | given | not mentioned |
| Read layout | paired-end HiSeq 2000 | TruSeq R1 / i7 index / R2 |

No contradiction found; upstream omits steps that do not change the library structure.
Whether the 2015 protocol changed cycle numbers or polymerases cannot be checked 🔴.

## 8. Open questions

- 🔴 Adapter concentration, ligation time/temperature; how the methylated adapters were
  obtained (bought pre-methylated vs. synthesized with 5-methyl-dC).
- 🔴 Primer identity and template for each of the two PCR rounds (2013) and whether the
  2015 protocol keeps two rounds.
- 🔴 Read lengths; the index set used (which 6-nt TruSeq indices).
- 🔴 Supplemental Table 5 nested-PCR primers (sperm MspI/HpaII validation).

## 9. How this note was made (tool evaluation)

`tools/get_sources.py` fetched the PMC full text, the upstream page and the two
*Nature Protocols* supplements; it did not fetch (or list as manual) the *Nature
Protocols* main text or the Genome Research supplement. `tools/scrape_primers.py` found
sequences **only** in the upstream page — correctly, since neither primary text fetched
contains any oligo; its `lib/` annotations (P5, TRUSEQ_READ1, INDEX1_PRIMER, P7_RC)
matched the computation above. The methods were read by hand.
