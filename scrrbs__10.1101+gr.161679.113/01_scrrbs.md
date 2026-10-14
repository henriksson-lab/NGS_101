# scRRBS — single-cell reduced-representation bisulfite sequencing

## Adapter-source audit

🟢 The defining paper calls the ligated oligos “Illumina standard premethylated indexed
adaptors,” and uses paired-end HiSeq 2000 sequencing. Illumina's adapter-sequence guide
publishes the matching TruSeq universal adapter, six-base single-index adapter, Read 1,
Index 1 and Read 2 primers. The library is single-indexed; it does not contain i5.

> **Evidence marking.** 🟢 verbatim from the source · 🟡 derived or inferred · 🔴 not
> published / not available. Relationships marked 🟡 *(computed)* were worked out with
> `lib/`; the model is `tools/scrrbs.py` on top of `lib/rrbs.py`.
>
> **The scg_lib_structs drawn page is not a source here.** It is inspiration only: no
> sequence, length or reaction in this note or in the model is taken from it, and
> agreement with it is not evidence. Every adapter and primer base drawn for scRRBS comes
> from `lib/illumina.py`, i.e. from *Illumina Adapter Sequences*, Document # 1000000002694
> — the vendor document for the "Illumina standard premethylated indexed adaptors" the
> paper says it ligated.

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
| `upstream_scRRBS.html.txt` | scg_lib_structs page | **inspiration only, not a source** — nothing is cited from or verified against it |
| `lib/illumina.py` (in-repo, from *Illumina Adapter Sequences* #1000000002694 v22) | canonical TruSeq universal / single-index adapter, Read 1, Read 2 and Index 1 primers | **every base the model draws** 🟢 (as Illumina publishes them) |
| `nprot.2015.039_…MOESM430_ESM.pdf.txt` | *Nat Protoc* Supplementary Figures 1–5 | nothing chemical (CpG coverage, methylation histograms) |
| `nprot.2015.039_…MOESM429_ESM.zip` → `zip/*.sh`, `SingleC_MetLevel.pl` | *Nat Protoc* analysis scripts | trimming / alignment parameters (§6) |

Not obtained — fetch by hand:

| What | URL | Why it matters |
|---|---|---|
| *Nature Protocols* main text (paywalled) | https://www.nature.com/articles/nprot.2015.039 | the only source that would name the library PCR primers, the index set, the reagent table and the 2015 cycle numbers |
| Genome Research supplement (not linked from PMC) | https://genome.cshlp.org/content/23/12/2126/suppl/DC1 | Supplemental Table 5 (nested-PCR primers for the MspI/HpaII sperm test), Supplemental Figs. 9–11 |

The Genome Research paper itself **gives no oligo sequence**: it says only "Illumina
standard premethylated indexed adaptors", ligated to dA-tailed MspI fragments with T4
ligase, and sequenced paired-end on a HiSeq 2000 🟢 (Methods, `PMC3847781.html.txt:1226`).
That phrase is a *pointer to a vendor document*, and the vendor document is authoritative
for the bases: this is the same provenance shape as scDNase-seq (paper → named standard
Illumina kit → Illumina's adapter document). What stays 🟡 is the identification — that
"standard premethylated indexed" means the TruSeq universal + 6-nt single-index pair —
and 🔴 the index set actually used and the premethylation sourcing.

---

## 1. What it is

Bulk RRBS (Meissner/Gu: MspI digest → end repair/A-tail → ligate methylated adapters →
bisulfite → PCR) shrunk to one cell. The change is not new chemistry but **removing every
purification before bisulfite**: lysis, lambda spike-in, MspI digestion, fill-in/A-tailing,
adapter ligation and the bisulfite reaction itself all happen **in the same tube** 🟢; the
first clean-up is the Zymo column after bisulfite, with 10 ng tRNA as carrier 🟢.

No cell barcode, no UMI. One cell = one library, identified by the **6-nt TruSeq index in
the ligated adapter** 🟡 (the paper only says "indexed"; 6 nt is the standard TruSeq
single-index width). No reverse
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

The paper names the adapters by vendor designation, not by sequence. The model therefore
imports them from `lib/illumina.py` and retypes nothing (`lib/rrbs.py`, which scRRBS,
scMT-seq and scTrio-seq share):

| Role in scRRBS | Model expression | Bases |
|---|---|---|
| universal (P5-side) adaptor strand | `illumina.TRUSEQ_P5_FULL` = `P5` + `TRUSEQ_READ1[4:]` | `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT` 🟢 |
| indexed (P7-side) adaptor strand | `illumina.INDEX1_PRIMER` + 6-nt index + `illumina.P7_RC` | `GATCGGAAGAGCACACGTCTGAACTCCAGTCAC`[index]`ATCTCGTATGCCGTCTTCTGCTTG` 🟢 fixed parts, index 🔴 |
| 12-bp fork stem | `illumina.STEM` | `GATCGGAAGAGC` 🟢 |
| Read 1 / Read 2 / Index 1 read primers | `seqprimers.TRUSEQ["R1"/"R2"/"I1"]` | the stock TruSeq primers 🟢 |
| cell index | 6-nt `N` placeholder | 🔴 — which TruSeq index set was used is not published |
| library PCR primers | graph inputs without a sequence state | 🔴 — the paper names neither primer (see §3) |

🟡 The *premethylation* of the adapter cytosines is a property the paper asserts
("premethylated") but no base in it changes: 5-methyl-dC is drawn as C.

🔴 The *Nature Protocols* 2015 step-by-step version of this protocol is the one source
that would name the PCR primers and the index set; it is paywalled and we do not hold it,
so nothing is asserted about primer names here (an earlier revision of this note carried
the names "QP1"/"QP2", which reached it only through the drawn page).

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
- Because the adaptor already contains full P5 and P7, primers no longer than
  `illumina.P5` / `illumina.P7` suffice to amplify, and no index is added by PCR 🟡
  *(computed)*; the primers the authors actually used are not named anywhere we hold 🔴.
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
   `5'-CGG…CCGA-3'`. 🟡 (mechanism, from the enzymes and dNTPs the paper lists) The
   filled-in C is from
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
   72 °C 1 min). Primers for either round are not named in the paper, and the model's
   graph therefore gives its two library-PCR primer inputs no sequence state and marks the
   step's geometry unavailable 🔴. Whether round 2 is on a dilution or clean-up of round 1
   is not said 🔴.
8. **Size-select** 200–500 bp on a 12 % native polyacrylamide TBE gel; Fragment Analyzer
   QC; qPCR quantification.

In the first PCR cycle only the P7-side primer can anneal: each converted strand carries
the universal adaptor at its 5' end and the P7 complement at its 3' end, so a P5 primer
finds a site only on the copy. Both original strands give the same structure. 🟡
*(computed from the adapter design)*

Side assay (not part of the library): single sperm lysed the same way, digested with MspI
or the methylation-sensitive **HpaII**, then nested PCR (30 + 30 cycles). The paper
assigns its nested-PCR primers to Supplemental Table 5; these are validation primers,
not scRRBS library-construction oligos. 🔴 (table not obtained)

## 4. Final library — 🟡 (assembled from the canonical adapter parts)

```
5'- P5 · TruSeq Read 1 · <insert: Y-GG … Y-T-G> · A · INDEX1_PRIMER seq (= Read 2') · <i7, 6 nt> · P7' -3'
```

- `P5 · TruSeq Read 1` = the universal adaptor (58 nt); `INDEX1_PRIMER · i7 · P7'` = the
  indexed adaptor (63 nt). Adapters add **121 nt** (+ the A-tail) to the insert.
- Insert: starts at the MspI site, `CGG` if that CpG was methylated, `TGG` if not
  (bisulfite); ends `…C/T` `T` `G`, the `T` being the artificial filled-in C; then the `A`
  tail. Non-CpG C are expected to read as T.
- Size-selected 200–500 bp total 🟢 → insert roughly 80–380 bp 🟡.

## 5. Sequencing

🟢 Paired-end on an **Illumina HiSeq 2000**; read length not stated in the Genome
Research paper 🔴.

🟡 *(computed; standard TruSeq read primers)*:

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

## 7. What the paper fixes and what it leaves open

| Point | Genome Research 2013 | Status in the model |
|---|---|---|
| Adapters | "Illumina standard premethylated indexed adaptors", no sequence 🟢 | canonical TruSeq parts from `lib/illumina.py`; the identification is 🟡 |
| Enzyme / fill-in / A-tail | MspI (C^CGG), Klenow exo⁻, dATP + dGTP + dCTP 🟢 | drawn as a digestion → end-repair → dA-tail chain |
| PCR | **two rounds** (PfuTurbo Cx 25 cyc, then Phusion 22 cyc) 🟢, primers unnamed 🔴 | one PCR step with sequence-less primer inputs and unavailable geometry |
| Lambda spike-in, tRNA carrier, PAGE size selection | given 🟢 | conditions, not construct structure |
| Read layout | paired-end HiSeq 2000 🟢; read lengths not stated 🔴 | Read 1 / Index 1 / Read 2 roles, index width 🟡 |

Whether the 2015 *Nature Protocols* version changed cycle numbers, polymerases or primer
names cannot be checked without that paper 🔴.

## 8. Open questions

- 🔴 Adapter concentration, ligation time/temperature; how the methylated adapters were
  obtained (bought pre-methylated vs. synthesized with 5-methyl-dC).
- 🔴 Primer identity and template for each of the two PCR rounds (2013) and whether the
  2015 protocol keeps two rounds.
- 🔴 Read lengths; the index set used (which 6-nt TruSeq indices).
- 🔴 Supplemental Table 5 nested-PCR primers (sperm MspI/HpaII validation).

## 9. How this note was made (tool evaluation)

`tools/get_sources.py` fetched the PMC full text and the two *Nature Protocols*
supplements (and the drawn page, which is not used); it did not fetch (or list as manual)
the *Nature Protocols* main text or the Genome Research supplement.
`tools/scrape_primers.py` found **no oligo at all** in the sources that count: neither
fetched primary text contains one. That is the finding — the adapter bases are available
only from Illumina's own adapter document, which `lib/illumina.py` carries. The methods
were read by hand.
