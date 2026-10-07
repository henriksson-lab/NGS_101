# Datasets — CRISPR-UMI (Schmierer)

**Schmierer B, Botla SK, Zhang J, Turunen M, Kivioja T, Taipale J.** "CRISPR/Cas9 screening
using unique molecular identifiers." *Mol Syst Biol* 2017;**13**(10):945.
doi:[10.15252/msb.20177834](https://doi.org/10.15252/msb.20177834) · PMC5658704 · CC BY

Chemistry: [`../../crispr-umi-schmierer__10.15252+msb.20177834/01_crispr-umi-schmierer.md`](../../crispr-umi-schmierer__10.15252+msb.20177834/01_crispr-umi-schmierer.md)

## Availability 🟢

> "Raw sequencing data: European Nucleotide Archive, PRJEB18436. Computer scripts: GitHub
> http://github.com/zhjilin/RSLC"

| | |
|---|---|
| study | **PRJEB18436** (secondary **ERP020364**), submission ERA770577 |
| submitter | Dept of Biosciences and Nutrition, Karolinska Institutet |
| released | 2017-10-24, status PUBLIC |
| code | `github.com/zhjilin/RSLC` (the pipeline; it does **not** ship the library file it needs) |
| **library** | **Dataset EV1** = `MSB-13-945-s002.csv`, ✅ **23,279 unique guides** |
| runs | 5 |
| ✅ total | ~6.29 × 10⁹ reads, ~127 GB gzipped |

## ✅ The guide library is published — Dataset EV1

🟢 Methods: "The guide library targets 2,325 genes and contains a total of **23,279 guides**
(Dataset EV1) … all human transcription factors, other genes of interest as well as
ribosomal proteins as positive controls and **101 non-targeting guides** as negative
controls. All sgRNA sequences used in this library were taken from a previously published,
genome-wide library" (Wang *et al* 2015).

✅ The file is `MSB-13-945-s002.csv`, 869 KB, and is already in exactly the `lib.csv` format
the authors' RSLC pipeline expects — `GuideID,GuideSequence,TargetGene`, **no header row**:

```
AATF_01,GTAGATAGGAAAGCCAGCAA,AATF
NONT1_01,GGTCCATGGGTGGAGTTACG,NONT1
```

✅ 23,332 rows, 23,279 unique 20-mers (so 53 sequences are shared between gene entries).

⚠ It is **not** on the RSLC GitHub repo, which ships only the Perl/R scripts and expects you
to supply `lib.csv` yourself. Fetch it with `gcbias/download/get_libraries.py schmierer`;
PMC's direct `/bin/` URL is behind a JavaScript interstitial, so the script uses the
EuropePMC `supplementaryFiles` REST endpoint instead.

## The experiment, as the paper describes it 🟢

> "we screened the human colorectal carcinoma cell line RKO for essential genes with an
> RSL-guide library targeting 2,325 genes with 10 guides per gene (Wang *et al*, 2015).
> Briefly, Cas9-expressing RKO cells were transduced with the lentiviral guide library, and
> samples were taken at Day 4 and Day 28 after transduction (control and treatment time
> points, respectively) … The experiments were run in duplicate and at far larger screen
> size … and sequencing depth (reads per guide) than previous screens"

| | |
|---|---|
| cell line | RKO (human colorectal carcinoma), ATCC; Cas9 via `pLenti-Cas9-sgHPRT1` |
| library | 2,325 genes × 10 guides ≈ **23,250 guides**, from Wang *et al* 2015 |
| | paper elsewhere says "> 23,000 guide-sets" 🟢 |
| timepoints | **Day 4** (control) and **Day 28** (treatment) |
| replicates | 2 |
| RSL | **6 nt**, cloned into the library |
| controls | ribosomal proteins (positive), non-targeting (negative), MYC called out separately |
| ✅ depth | ~2.7 × 10⁵ reads per guide summed over runs — very deep |

🟢 Library coverage in cells: "Both time points together covered 93% of input RSL-guides
(**78 million unique sequences** in the cell population)". D4∩D28 overlap ≈ two-thirds.

🟢 Analyses defined in the paper: **TCA** (total count), **LDA** (lineage dropout — fraction
of RSL-guides lost D4→D28), **IRA** (internal replicate — RSL-guides binned into **64**
internal replicates). Subsampling to ¼ and ¹⁄₁₆ was used to simulate smaller screens.

## The runs

🟡 Sample titles are ENA's; the paper does not print a run table, so the mapping below is the
obvious reading of the titles against the design. `NTU` is 🔴 unexplained in the paper.

| run | sample title | author's filename | ✅ reads | ✅ mean len | gz |
|---|---|---|---|---|---|
| ERR2114695 | `RKO_NTU_INPUT` | `input.fq.gz` | 2,503,174,977 | 23.6 | 49 GB |
| ERR2114696 | `RKO_NTU_Replicate1_D4` | `ntu_r1_d4.fq.gz` | 660,790,803 | 27.0 | 14 GB |
| ERR2114697 | `RKO_NTU_Replicate1_D28` | `ntu_r1_d28.fq.gz` | 1,189,493,924 | 22.7 | 24 GB |
| ERR2114698 | `RKO_NTU_Replicate2_D4` | `ntu_r2_d4.fq.gz` | 794,719,413 | 27.0 | 17 GB |
| ERR2114699 | `RKO_NTU_Replicate2_D28` | `ntu_r2_d28.fq.gz` | 1,145,941,439 | 22.9 | 23 GB |

All runs: `library_strategy=OTHER`, `library_layout=SINGLE`, Illumina HiSeq 4000.

## ✅ Where the RSL actually is — read this before downloading

The paper's readout is `20 + 6 + 6` cycles, with **i7 reading the RSL** and i5 the sample
index. That means the RSL is in an **index read**, and index reads are handled differently by
the two files ENA serves for each run:

| file | header | RSL present? |
|---|---|---|
| ENA-generated `ERR*.fastq.gz` | `@ERR2114696.1 K00110:106:…/1` | ✅ **no** — indexes stripped |
| author-submitted `ntu_*.fq.gz` | `@K00110:106:…:1261 1:N:0:CCACTT+ATCCCG` | ✅ **yes** — in the index field |

> ⚠ **Use the `submitted_ftp` files, not the `fastq_ftp` files.** The generated FASTQ renames
> every record and discards the index pair, which is the entire UMI. Nothing in the ENA
> metadata warns you about this.

Paths follow `ftp.sra.ebi.ac.uk/vol1/run/ERR211/<run>/<filename>`, e.g.
`…/ERR2114696/ntu_r1_d4.fq.gz`.

✅ Verified on the first ~123,000 reads of `ntu_r1_d4.fq.gz`:

```
header index field = i7 + i5     e.g. 1:N:0:CCACTT+ATCCCG
  i7  4,056 distinct values  -> the 6-nt RSL   (4^6 = 4,096 possible)
  i5     18 distinct values  -> the sample index
```

The assignment is confirmed empirically, not just taken from the text: i7 varies read to read
as a random label must, while i5 is dominated by a couple of fixed values.

✅ Read 1 is the 20-nt spacer plus a few low-quality run-off cycles (mean length 22.7–27.0,
with the tail inconsistent between reads and quality strings like `FF-A--`). Treat anything
past base 20 as junk.

## ⚠ Two problems to settle before trusting this data

**1. The RSL space is saturated.** See [`README.md`](README.md) — ✅ 81.9 % of the 95.2 M
possible (guide, RSL) pairs are occupied, leaving only ~68 % of lineages within a guide
distinguishable. The paper additionally filters low-read-count RSLs 🟢, so the effective
denominator is not even the full space. **Distinct-RSL counts are not molecule counts here.**

✅ The observed i7 distribution is also not flat: the most common values are C-rich
(`CCCCCC` 826, `ACCCCC` 596, `CCACCC` 407, `GCCCCC` 369, `TCCCCC` 355 in that 123 k sample).
🔴 Whether that is amplification bias, index-read base-calling behaviour, or skew in the
cloned RSL pool is **not resolved** — and it matters a great deal, because a composition skew
in the label itself would contaminate any GC-bias estimate built on those labels. Resolve it
before using this data for the reviewers' question.

**2. The submitted file may not be one sample.** ✅ In `ntu_r1_d4.fq.gz` the i5 field shows
`ATCACG` (~70 %) and `ATCCCG` (~29 %) plus 16 rare variants. Those two differ by one base, so
this is either two pooled samples or systematic miscalling of one index — 🔴 undetermined.
Check before assuming a run equals a sample.

## What it is good for

🟡 Despite the above, the **`RKO_NTU_INPUT`** run is the single most interesting sample in
either paper for the GC question: it is plasmid/input material, so a GC trend in it cannot be
explained by differential growth, and at 2.5 × 10⁹ reads it is extremely deep. The RSL
saturation argument applies to using RSLs as a *molecule counter*; it does not prevent using
read counts per guide against guide GC%, with the input as the no-selection baseline.
