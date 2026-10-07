# Datasets — CRISPR-UMI (Michlits)

**Michlits G, Hubmann M, Wu S-H, Vainorius G, Budusan E, Zhuk S, *et al*.** "CRISPR-UMI:
single-cell lineage tracing of pooled CRISPR–Cas9 screens." *Nat Methods*
2017;**14**(12):1191–1197. doi:[10.1038/nmeth.4466](https://doi.org/10.1038/nmeth.4466)

Chemistry: [`../../crispr-umi-schmierer__10.15252+msb.20177834/03_crispr-umi-michlits.md`](../../crispr-umi-schmierer__10.15252+msb.20177834/03_crispr-umi-michlits.md)

## Availability 🟢

> "All genomic data are available via the NCBI Sequence Read Archive under accession code
> PRJNA383356."

| | |
|---|---|
| study | **PRJNA383356** |
| runs | 4 |
| ✅ original per-lane BAMs | 9 |
| internal project name | `CRISPR_DigDeep` (appears throughout the SRA titles) 🟢 |

## The library 🟢

| | |
|---|---|
| barcode | **10 random nt**, introduced by parallel cloning into a retroviral backbone |
| genes | 6,560 mouse genes — all genes encoding nuclear proteins |
| guides | 4 sgRNAs/gene + **112 non-targeting controls**; full library **26,514 guides** |
| | ⚠ ✅ 6,560 × 4 + 112 = **26,352**, which is 162 short of the stated 26,514 (0.6 %). 🟡 Subpools were cloned separately, so the extra are most likely additional controls or a few genes with 5 guides; 🔴 not stated. Use the library file, not the arithmetic. |
| cloning coverage | **954–8,776** independent cloning events per sgRNA |
| complexity | **83.5 million** UMIs — "exceeds the number of clones assayed in a screen" |
| representation | 10th-to-90th-percentile spread in the library = **fourfold** |
| design | Doench score, off-target, position in transcript, exon structure, Pfam domains |
| pilot library | 365 genes / **1,437 guides** |

🟢 Enrichment before PCR: PacI digestion plus size selection enriches the cassette
**10³–10⁴-fold**, and the stated purpose is directly relevant here — "This enrichment
minimized the number of PCR cycles required for amplification, **thereby reducing PCR
amplification biases**." The 589-bp PacI band is gel-excised (control amplicon 582 bp on a
7,710-bp PacI fragment).

## The runs

🟡 Titles below are the SRA `TITLE` fields, which are more informative than the `sample_title`
and match the paper's three screens plus the library check.

| run | SRA title 🟢 | ✅ reads | ✅ BAMs |
|---|---|---|---|
| SRR5559298 | `CRISPR_DigDeep_KO_Library: Plasmid Library Sequencing` | 625,912,869 | 2 |
| SRR5559297 | `…KO screen: Ms embryonic stem cells, pilot screen hypersensitivity to Etoposide` | 292,762,401 | 2 |
| SRR5559299 | `…KO_screen: Ms embryonic stem cells, Hypersensitivity to Etoposide, Negative Selection` | 563,804,278 | 4 |
| SRR5559300 | `…KO_screen: Ms induced pluripotent stem cells, Roadblocks of reprogramming MEF to iPSC, Positive Selection` | 73,462,865 | 1 |

All runs: single-end, exactly 50 nt, Illumina HiSeq 2500.

### Design of each screen, from the paper 🟢

- **Negative selection, etoposide hypersensitivity** — mESC with Dox-inducible Cas9, sgRNAs
  delivered retrovirally. **4,000 independently infected cells per sgRNA**. Only **2 samples**
  (control and etoposide-treated); the screen was "carried out once". The cell number was
  chosen from "the expected number of reads obtained from a SR50 HiSeq run (about 200 million)
  and the size of the library (about 26500 guides) and the number of samples (2)".
- **Limiting-dilution arms** in the pilot: (i) 20 clones of 200 cells, (ii) 64 clones of
  64 cells, (iii) 200 clones of 20 cells, (iv) no dilution → 4,000 mostly independent cells.
- **Positive selection, MEF → iPSC reprogramming** — **4 replicates from 4 embryos**;
  Oct4-GFP⁺ cells sorted per replicate. MEFs from E13.5.

🟢 Barcode error handling: "adjacency-based barcode collapse, given an edit distance of 1 and
a threefold read-count difference". Thresholds: clones with ≥3 total reads, guides with ≥5
clones present. Reads assigned with **Bowtie**.

## ✅ The library, primers and index map are published

The article is not open access and has no PMC record, but its supplementary objects are
served without a paywall from Springer static-content. `gcbias/download/get_libraries.py`
fetches all three:

| table | file | contents |
|---|---|---|
| Suppl. Table 2 | `MOESM3_ESM.xlsx` (2.0 MB) | the **guide library** — `group`, `20nt guide`, `oligo name`, `77nt oligo`. ✅ 28,565 rows, **26,486 unique** 20-mers (paper text says 26,514) |
| Suppl. Table 3 | `MOESM4_ESM.xlsx` | PCR and sequencing primers (`FW_G_CrSc…`, `RV_G_CrSc…`) |
| Suppl. Table 4 | `MOESM5_ESM.xlsx` | **sample ↔ 6-bp experimental index** map — `Sample number`, `Barcode`, `Treatment`, `Screen-Format`, plus explicit `Bamfile:` entries naming the original BAM per sample |

🟡 Supplementary Table 4 is the missing half of the BAM route below: it says which 6-bp index
belongs to which arm (`ctrl` vs `etopos`, the limiting-dilution formats, `none - plasmid
lib`), which is exactly what cannot be recovered from the public FASTQ. It is only useful if
the original BAMs are obtained.

## ⚠ The UMI is NOT in the public FASTQ

This is the decisive fact for reuse. 🟢 The SRA `DESIGN_DESCRIPTION` states the layout:

> "SR50, DUAL INDEXING, Read 1: custom primer (first 20bp are sgRNA Sequence) **index read1:
> UMI (first 10bp random nucleotides used as unique molecular identifiers) index read2: 6bp
> sequences for demultiplexing**"

✅ But the SRA run itself reports `nreads="1"`, one 50-nt read per spot. Confirmed by reading
the file: records are `@SRR5559298.1 1/1` with no index field, and bases 21–50 are
predominantly `N` with a fixed pattern of stray calls — i.e. run-off cycles past the 20-nt
guide, not data.

```
@SRR5559298.1 1/1
NGAGAACTCCAACTTCCCCGNNNTNAGANNTNNNNNNGNNNNAGNNTAGN
+
#<<BB//<F//B//BF/F<F###B#<<F##<######<####<<##<<<#
 |<--- 20-nt sgRNA --->|<------ run-off, not sequence ------>|
```

Consequences 🟡:
1. **No UMI** → no clonal/lineage analysis, no deduplication.
2. **No experimental index** → the pooled samples inside a run *cannot be separated*. Most
   importantly, SRR5559299 contains both the control and the etoposide-treated arm with no way
   to tell them apart from the FASTQ.
3. There is also **no author-submitted FASTQ** (`submitted_ftp` is empty) — unlike Schmierer,
   where the index survives in the submitted file's header.

### ✅ The index reads probably do exist — in the original BAMs

SRA retains the submitters' original per-lane BAMs, which for this kind of core-facility
output normally carry the index reads as tags. They are listed as `supertype="Original"`:

| run | original BAM | size |
|---|---|---|
| SRR5559298 | `Plasmid_NDS_C89KRANXX_6_20160127B_20160128.bam` | 7.36 GB |
| SRR5559298 | `Plasmid_NDS_C9PLEANXX_1_20161111B_20161114.bam` | 12.58 GB |
| SRR5559297 | `Pilot_screen_C9B8KANXX_7_20160406B_20160408.bam` | 2.07 GB |
| SRR5559297 | `Pilot_screen_C8LABANXX_2_20160414B_20160416.bam` | 7.08 GB |
| SRR5559299 | `30k_Etoposide_CAF25ANXX_1_07-02-2017.bam` | 6.89 GB |
| SRR5559299 | `30k_Etoposide_CAF25ANXX_2_07-02-2017.bam` | 3.93 GB |
| SRR5559299 | `30k_Etoposide_CAR6CANXX_7_20170310B_20170311.bam` | 4.84 GB |
| SRR5559299 | `30k_Etoposide_CAR6CANXX_8_20170310B_20170311.bam` | 4.23 GB |
| SRR5559300 | `30k_iPSC_H5F3TBCXX_1_20160927B_20160928.bam` | 2.90 GB |

> 🔴 **Not anonymously downloadable.** They are flagged `access_type="Use Cloud Data Delivery"`
> with `free_egress="-"`; ✅ a direct HTTPS range request to the `sra-pub-src-*` bucket returns
> **403**, and the NCBI resolver offers only the ETL-processed run (which is the UMI-less
> FASTQ). Getting them means NCBI Cloud Data Delivery, i.e. a requester-pays S3 transfer of
> ~52 GB — or asking the authors directly, which is likely faster.
>
> 🟡 It is also **unverified** that the index reads survived into those BAMs. Confirm on the
> smallest one (`30k_iPSC…`, 2.90 GB) before paying to move the rest.

## ⚠ A contradiction between figure and methods

🟢 Figure 2c labels the construct `A5 – Index – U6 – sgRNA – i7 – NNNNNN – A7` with
"Index 1: experimental index" and "Index 2: barcode".

🟢 The Methods say the opposite: "it is necessary to obtain at least **10 bp for index 1
(barcode)** and **6 bp for index 2 (experimental index)**" — and the SRA
`DESIGN_DESCRIPTION` agrees with the Methods (index read 1 = UMI, index read 2 = 6 bp
demultiplexing).

🟡 **The Methods and SRA are self-consistent and the figure is the odd one out**, so take
index 1 = the 10-nt UMI. Figure 2c also draws the barcode as `NNNNNN` (6 N) where the text
says 10 — consistent with the figure being schematic rather than exact. Worth knowing before
writing a parser against the figure.

🟢 Custom read primer (matches `crispr-umi-schmierer__10.15252+msb.20177834/tools/michlits.py`):
`CGATTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCG`

## What it is good for

🟡 If the original BAMs do carry the index reads, **this is the better of the two datasets**
for the GC question: ✅ the (guide, UMI) space is only **0.300 %** occupied
(arithmetic and the comparison against Schmierer in [`README.md`](README.md)), so ~99.9 % of
lineages are separable and distinct-UMI counts behave like real molecule counts. The
`Plasmid Library Sequencing` run is then the ideal no-biology baseline, and the PacI
enrichment means it was amplified with deliberately few cycles — exactly the contrast
CRISPR-MIP is arguing about. Until the BAMs are checked, none of that is usable.
