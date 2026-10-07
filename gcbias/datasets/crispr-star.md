# Datasets — CRISPR-StAR (Michlits / Elling lab, 2025)

**Ueberschär E, Sommer A, … Elling U.** "CRISPR-StAR enables high-resolution genetic
screening in complex in vivo models." *Nat Biotechnol* 2025.
doi:[10.1038/s41587-024-02512-9](https://doi.org/10.1038/s41587-024-02512-9) ·
PMC12611787 · **open access**

🟡 Found while looking for a dataset with a *usable* lineage UMI. It is the follow-up to
[CRISPR-UMI (Michlits)](crispr-umi-michlits.md) and cites it for the UMI design, so it
carries the same **10-nt** barcode — but, unlike the 2017 paper, the barcode is in a **read**
rather than a stripped index.

## Why it matters here

The two 2017 CRISPR-UMI papers fail for opposite reasons ([README](README.md)): Schmierer's
RSL is readable but 82 % saturated, Michlits' is sparse but was never deposited. ✅ This one
appears to have both properties at once:

| | Schmierer | Michlits 2017 | **CRISPR-StAR** |
|---|---|---|---|
| UMI length | 6 nt | 10 nt | **10 nt** |
| UMI in public data | yes (FASTQ header) | **no** | ✅ **yes** — and already as counts |
| saturation | 82 % | 0.3 % | 🟡 expected ~sparse |
| platform | HiSeq 4000 | HiSeq 2500 | **NextSeq 2000** |
| internal control | no | no | ✅ **yes** — active/inactive sgRNA in the same clone |

🟢 The platform difference is independently useful: the GC deficit measured so far is on
HiSeq 4000 only, and NextSeq 2000 uses two-colour chemistry, so agreement across them would
argue against a platform artefact.

## Availability 🟢

> "All primary and processed data of the genome-wide CRISPR–Cas9 screen are uploaded to the
> Gene Expression Omnibus database under series **GSE262309**. As the screen was performed
> in two batches, batch 1 is available under **GSE262307**, and batch 2 under **GSE262308**."

🟢 Read structure, from the GEO `data_processing` field:

> "Read 1 (R1) contained the sgRNA sequence and to confirm an active or inactive sgRNA, we
> assessed for each Read 1 sequence the last 6 bases. Reads that contained 'TTTT' in the end
> were assigned as inactive sgRNA and 'CAGC' containing reads were marked as active sgRNAs.
> **Read 2 (R2) contained the UMI barcode sequence.**"

Instrument: NextSeq 2000. `library_strategy = OTHER`.

### ✅ The per-sample files are already sgRNA-UMI count tables

No FASTQ parsing needed. Each GSM carries two supplementary files, `*_active_MD_noShadows.txt.gz`
and `*_inactive_MD_noShadows.txt.gz`:

```
guide                index                   UMI          reads
1110032F04Rik_1      in_vitro_Blue_active    GAATGAATTA      38
1110032F04Rik_1      in_vitro_Blue_active    AACGGATAAC      25
```

✅ Verified by download: `GSM8163104_Lane1_in_vitro_Blue_active_MD_noShadows.txt.gz`,
13.2 MB, **2,162,437 rows**, 10-nt UMIs.

### 🔴 Do NOT use these tables for a GC measurement — use the raw FASTQ

`MD_noShadows` is not a neutral tidy-up. 🟢 Methods, on how "shadows" are removed:

> "sorted the dataset in descending order based on the total number of reads associated with
> each UMI … iteration through the UMIs with the highest number of reads until the total
> number of reads associated with a UMI drops below the predefined threshold of 100,000 …
> A ratio was calculated for each occurrence by dividing its number of reads by the maximum
> number of reads among all occurrences of the same UMI. These ratios were then used to
> **filter out instances where the ratio was less than or equal to 0.001**"

⚠ That is a **read-count-dependent filter**, and read count is exactly the quantity under
measurement. A guide that amplifies poorly has fewer reads per (guide, UMI) occurrence, so
it is more likely to fall under the 0.001 ratio and be deleted — which would inflate an
apparent GC effect by construction. The filter is sensible for its own purpose (removing
UMI–guide chimeras and index hopping) and wrong for ours.

### ✅ The raw FASTQ carries the UMI — this is what to use

✅ Verified on `SRX24034750` / `SRR28430769` (GSM8163104):

| | |
|---|---|
| reads per spot | **2** — unlike Michlits 2017, where it was 1 |
| read 1 | 75 nt — the sgRNA, plus the last 6 bases encoding active (`CAGC`) / inactive (`TTTT`) |
| read 2 | **10 nt — the UMI** |
| spots | 28,487,791 |
| original files | `Lane1_in_vitro_Blue.R1.fastq.gz` (1.05 GB), `.R2.fastq.gz` (0.40 GB) |

🟢 So the UMI is in a genuine read in the public archive, and the whole analysis can be run
from raw reads with our own error-collapse rather than theirs. The processed tables remain
useful as a cross-check — if our numbers and theirs disagree, the filter above is the first
thing to look at.

## 🟡 The guide library — mapping names to sequences

The count tables name guides `<gene>_<n>`, not by sequence, so GC% needs the library.

| | file | contents |
|---|---|---|
| Suppl. Table (MOESM3) | `41587_2024_2512_MOESM3_ESM.xlsx`, 7.8 MB | **mouse** library — `ensembl_gene_id`, `external_gene_name`, `description`, `line.sublist`, `GeCKO.cloning.FW_oligo`, `GeCKO.cloning.RV_oligo`, **`sgRNA`**; 4 sheets (`A1.drugged`, `A2.surface`, `A3.metabol_catabol`, `A4.drug.rest`), ✅ 105,637 20-mers |
| Suppl. Table (MOESM4) | `..._MOESM4_ESM.xlsx`, 7.3 MB | **human** library (ENSG ids), 100,370 20-mers |

🟡 **The join is plausible but unvalidated.** Sheet `A1` holds 15,584 guides over 3,115 genes
at **5 guides per gene** (3,112 of them), and the count-table suffixes run 1–5 with ~2,255
guides each — consistent with `<gene>_<n>` meaning the n-th guide of that gene in library
order. 🔴 That ordering assumption has **not** been checked against anything independent, and
getting it wrong would scramble GC assignment without any obvious symptom. Validate it
before use — e.g. against a gene whose guides differ enough in GC to be identifiable, or by
re-deriving counts from the raw FASTQ for a handful of guides.

⚠ Also note the library is **GeCKO-derived** cloning oligos, so the amplicon and primers
differ from both 2017 papers.

## What to do with it

🟡 In order of value:

1. **Repeat the GC measurement on a sparse UMI.** Schmierer's 6-nt RSL forced a collision
   correction that the data then rejected; a 10-nt UMI at genome-wide scale should make
   distinct-UMI counts behave like real molecule counts, so `reads per molecule` becomes
   nearly assumption-free.
2. **Cross-platform check.** NextSeq 2000 vs HiSeq 4000.
3. **The internal control is the real prize.** Active and inactive sgRNAs sit in the *same
   clone*, distinguished only by the last 6 bases of read 1 (`CAGC` vs `TTTT`). That is a
   paired comparison with identical template, identical clone, identical amplification —
   closer to a controlled amplification contrast than anything in the other datasets.
