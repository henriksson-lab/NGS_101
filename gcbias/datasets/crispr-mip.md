# Datasets — CRISPR-MIP (ours) and the external sets used for GC bias

**Preprint:** bioRxiv [2024.03.28.587082](https://www.biorxiv.org/content/10.1101/2024.03.28.587082v1.full)
Chemistry: [`../../crispr-mip__10.1101+2024.03.28.587082/01_crispr-mip.md`](../../crispr-mip__10.1101+2024.03.28.587082/01_crispr-mip.md)

## Our own data 🔴

🟢 Data availability, verbatim:

> "The raw sequencing data will be deposited on ArrayExpress. The pipeline for analyzing
> CRISPR-MIP sequencing data is available on GitHub,
> https://github.com/henriksson-lab/crispr-mip"

| | |
|---|---|
| raw data | 🔴 **no accession yet** — "will be deposited" |
| pipeline | `github.com/henriksson-lab/crispr-mip` 🟢 |

> ⚠ Action item independent of the GC-bias question: reviewers who ask for more analysis will
> also expect the accession to exist. "Will be deposited" does not survive a revision round.

## The GC-bias analysis already in the preprint 🟢

Worth having written down, because the reviewers are asking to extend *this*, and any new
analysis should not silently duplicate it.

> "**Analysis of GC% bias.** We searched for other CRISPR screens that might be GC-biased
> primarily using Pubmed and BioGRID ORCS. The processed Sintov2022 sgRNA counts were
> downloaded from GEO #GSM6008413. For the DepMap CRISPR screen, raw counts were obtained
> from https://score.depmap.sanger.ac.uk/downloads, Release 1 (5th April 2019),
> raw_sgrnas_counts.zip. The GC% bias was analyzed in R. Jitter has been added to the GC% in
> the plots, but models are fitted to the original values. For Sintov2022, the displayed
> linear model was gen…"

| dataset | accession / source | form | 🟡 note |
|---|---|---|---|
| **Sintov 2022** | GEO **GSM6008413** | *processed* sgRNA counts | a single GSM — one sample, already collapsed |
| **Sanger SCORE / DepMap** | `score.depmap.sanger.ac.uk/downloads`, Release 1 (2019-04-05), `raw_sgrnas_counts.zip` | **raw** sgRNA counts | large, many cell lines |
| search tools | PubMed, **BioGRID ORCS** | — | ORCS is the right index for finding more screens |

🟢 Method detail already fixed in the preprint: the GC% trend is removed by fitting a
polynomial in GC% to the sgRNA log fold changes and subtracting it; per-gene fold change is
the fold change of the corresponding sgRNAs.

## 🟡 Where the gap is, and why the UMI papers were the natural place to look

Both sets above are **count tables keyed by sgRNA**. They can show that abundance correlates
with guide GC%, but they cannot separate the two mechanisms that produce that correlation:

| mechanism | needs |
|---|---|
| **amplification/PCR bias** — GC-rich guides amplify differently | a molecule counter independent of read count |
| **biology** — GC correlates with gene class, essentiality, expression | a no-selection baseline |

A lineage UMI is attractive because it supplies the first: reads-per-molecule is a direct
amplification readout, and the ratio of reads to molecules against GC% is the measurement the
reviewers are implicitly asking for. That is why the two CRISPR-UMI papers were checked
first — and [`README.md`](README.md) records why neither delivers it off the shelf.

🟡 What is still achievable without a working UMI, in increasing order of effort:

1. **Input-sample GC regression.** Schmierer's `RKO_NTU_INPUT` (2.5 × 10⁹ reads, plasmid
   input) and Michlits' `Plasmid Library Sequencing` run are both essentially biology-free.
   A GC trend there is amplification, not selection. This needs no UMI at all and is the
   cheapest new evidence available. See the per-paper documents.
2. **Cycle-number contrast.** Michlits' PacI enrichment exists specifically to cut cycle
   count 🟢; if a low-cycle and a high-cycle library of the same material can be compared, the
   GC slope difference is the bias directly.
3. **Schmierer RSL, used carefully.** Not as a dedup counter — the space is 81.9 % occupied —
   but conditioned on (guide, RSL) pairs with a fixed, low read count, where saturation bites
   least. ⚠ Must be done with the saturation model explicit, and the C-rich skew in the
   observed RSL distribution resolved first.
4. **Michlits original BAMs**, if the index reads survived deposition and the ~52 GB
   requester-pays transfer is worth it. This is the only route to a genuinely sparse UMI.

✅ Item 1 is now implemented and run — see [`..`](../README.md). On the
Schmierer plasmid input the GC effect is small (≈0.7 % fewer reads per molecule per 10 % GC)
but strongly non-linear, flat through the mid-range and concentrated above ~65 % GC. A large
part of the apparent effect before UMI error-collapse was sequencing error rather than
amplification.

## ✅ Our own data already contains a DNA-input series — and the tooling to analyse it

🟢 From `script_slurm/samplemeta.csv` in `github.com/henriksson-lab/crispr-mip` (66 samples,
columns `protocol / seq_platform / cells / ug_dna / library / samplename / treatment / time /
replicate`). The design, by library and input:

| library | sample | template | **µg DNA** | cells | padlock | PCR |
|---|---|---|---|---|---|---|
| Brunello-kinome | `plasmid_lib` | **plasmid** | **0.2** | — | 3 | 4 |
| Brunello-kinome | `t0` | gDNA | **1** | — | 3 | 3 |
| Brunello-kinome | `normoxia_3d`, `normoxia_14d`, `hypoxia_3d` | gDNA | **10.1** | 1.5 × 10⁶ | 3+1+3 | 6+6+3 |
| Brunello (full) | `digested_gDNA` / `non-digested_gDNA` | gDNA | **7** | 1 × 10⁶ | 2 / 9 | — / 5 |
| Brunello (full) | `plasmid_lib` | **plasmid** | **0.2** | — | — | 4 |
| GFPg1 (1 guide) | digested / non-digested | gDNA | 1 | 1.5 × 10⁵ | 1 / 1 | — |

🟢 Guide sequences for GC% are in the same repo: `script_slurm/lib/kinome.csv` (6,304 guides,
column `sgRNA Target Sequence`) and `lib/brunello.csv`.

### ✅ Yes, there is a way to do this — the pipeline already exists

`subsamp_mip.py` emits a **`(grna, umi)`** table, which is exactly the shape
`gcbias/tools/extract.py` produces and that `rarefy.py` / `stratified.py` / `gcbias.py`
consume. So the chain is:

1. per-sample `(guide, UMI, reads)` table → the `*.guide_umi.tsv` format;
2. GC% per guide from `kinome.csv`;
3. ⚠ **`rarefy.py --target-rpu` to equalise reads-per-UMI across the samples being compared**;
4. `stratified.py` for the model-free GC statistic, `paired.py` for a per-guide contrast.

🟡 One parameter to match: their UMI clustering uses `UMIClusterer(threshold=2)` (edit
distance **2**), where our Schmierer work used adjacency collapse at distance **1**. Collapse
depth changed the measured GC effect materially on Schmierer, so this must be held constant
across any comparison rather than inherited from whichever script ran.

### 🔴 The trap: DNA input is the *worst possible* variable to compare naively

This is the one comparison where the compression problem in [`../README.md`](../README.md)
bites hardest, because **input amount sets molecules per guide directly**:

| sample | template | input | ✅ molecules of each guide per reaction |
|---|---|---|---|
| kinome `plasmid_lib` | plasmid | 0.2 µg | **≈ 2,260,000** |
| kinome `t0` | gDNA | 1 µg | **≈ 24** |
| kinome `3d` / `14d` | gDNA | 10.1 µg | **≈ 243** |
| Brunello `gDNA` | gDNA | 7 µg | **≈ 14** |

⚠ At a fixed sequencing depth, more molecules per guide means **fewer reads per molecule**,
and the compression factor on the GC statistic follows:

| reads per UMI | compression |
|---|---|
| 1.2 | **0.18** |
| 2.0 | 0.59 |
| 3.0 | 0.82 |
| 5.0 | 0.97 |

🔴 So a 10× input change at fixed depth moves compression by **2–3× on its own**, and a naive
analysis would report *"more DNA input → less GC bias"* as a **pure artefact**. The plasmid
samples are worse again: at ~2.3 × 10⁶ molecules per guide they sit roughly **10,000×** away
from the gDNA samples and are not on the same scale at all.

> **Rule:** rarefy every sample to a common reads-per-UMI *before* comparing, and report the
> operating point. Without that, this analysis measures the statistic, not the chemistry.

### 🟡 Which comparisons are actually matched

Taking the user's caution seriously — the comparison must hold everything but input constant:

| contrast | input change | what else moves | verdict |
|---|---|---|---|
| kinome `t0` (1 µg) **vs** `normoxia_3d` (10.1 µg) | **10×** | 3 days of culture, minimal selection | 🟢 **the usable one** |
| kinome `plasmid_lib` (0.2 µg) **vs** `t0` (1 µg) | 5× | plasmid vs gDNA template | 🟡 confounded by template type |
| Brunello `digested` **vs** `non-digested` (both 7 µg) | none | the digest itself | ✅ clean test of *accessibility*, not input |
| any kinome gDNA **vs** `plasmid_lib` | 50× | template type **and** ~10⁴× molecules/guide | 🔴 not interpretable |

### 🔴 There is NO sequencing data for an input titration — but there is one for probe

✅ Checked every row of `samplemeta.csv`: **`ug_dna` never varies within a comparable group.**

| group | input | varies? |
|---|---|---|
| GFPg1 digested vs non-digested | 1 µg both | no |
| Brunello gDNA — all 11 padlock + 5 PCR | **7 µg throughout** | no |
| Brunello-kinome screen samples | 10.1 µg throughout | no |
| …`t0` | 1 µg | — |
| …`plasmid_lib` | 0.2 µg | — |

🔴 The 0.2 / 1 / 7 / 10.1 µg spread tracks **what each sample is** (plasmid library vs t0 gDNA
vs screen gDNA), not a deliberate titration. So the DNA-input series exists **only as qPCR**
(Supplementary Fig 1a) and **was never sequenced**. The input-amount GC analysis sketched above
is therefore not runnable on existing data — it needs new libraries.

### ✅ What *was* sequenced is a 100× PROBE-concentration series

🟢 Brunello, padlock, non-digested gDNA, **7 µg and 10⁶ cells throughout**, same platform,
3 replicates each:

| probe | samples | ✅ probe molecules (20 µl) | ✅ probe : target |
|---|---|---|---|
| **0.2 µM** | `non-digested_gDNA_0.2_rep1-3` | 2.4 × 10¹² | 2,271,000 × |
| **0.02 µM** | `non-digested_gDNA_0.02_rep1-3` | 2.4 × 10¹¹ | 227,000 × |
| **0.002 µM** | `non-digested_gDNA_0.002_rep1-3` | 2.4 × 10¹⁰ | 22,700 × |

matching the preprint's stated range ("0.2 - 0.002 µM probe"). Plus `digested_gDNA_rep1-2` at
the same 7 µg, giving the digest contrast at fixed input.

🟡 **This is arguably a better test of the hypothesis than the input series would have been.**
Probe stays in vast excess even at the lowest concentration (≈22,700× over target sites), so
this is **not** stoichiometric limitation — it is a **hybridisation-kinetics** series. Lower
probe concentration means slower annealing during the overnight 60 °C ramp, and what sets that
rate is arm Tm and duplex stability, i.e. **sequence composition**. If capture is GC-biased,
the bias should *grow* as probe is diluted and the reaction moves further from completion.

✅ And capture is the step that is unique to CRISPR-MIP and where our chemistry is most
exposed — the Phusion fill-in at 60 °C, below its optimum, discussed in
[`polymerases.md`](polymerases.md). A GC trend that strengthens from 0.2 → 0.002 µM would
localise the effect to **hybridisation/capture** rather than to the PCR that follows.

> **So the runnable experiment on existing sequencing data is the probe titration, not the
> input titration** — three concentrations × three replicates, everything else held constant.
> Same pipeline: rarefy to a common reads-per-UMI, then `stratified.py`.

## ✅ RUN — the data is deposited, and both questions are answered

🟢 The raw data **is** public: **ArrayExpress `E-MTAB-14379`**, "CRISPR-MIP replaces PCR and
reveals GC and oversampling bias in pooled CRISPR screens". Deposited after the preprint, which
still says "will be deposited" — ⚠ that sentence should be updated to the accession.

Fetched with `gcbias/download/get_crisprmip.py` (4M read pairs per sample), extracted with
`gcbias/tools/extract_mip.py`. Runs `ERR13510408–418`.

### ✅ Result 1 — probe concentration does NOT change the GC effect

High-GC (>0.70) deficit, all nine samples rarefied to a matched **7.06 reads/UMI**:

| probe µM | rep1 | rep2 | rep3 | mean |
|---|---|---|---|---|
| 0.002 | −0.031 | −0.083 | −0.025 | **−0.046** |
| 0.02 | −0.059 | −0.069 | −0.018 | **−0.049** |
| 0.2 | −0.089 | −0.042 | −0.042 | **−0.057** |

✅ Pooled **−0.051 ± 0.024**. Between-group range **0.011**, within-group sd **0.024** — the
replicate scatter is more than double the spread across a **100×** probe range. 🔴 **No
resolvable trend.**

✅ Re-run aggregate-then-log (see the counting rule in [`../README.md`](../README.md)), the
means become −0.030 / −0.046 / −0.046, range 0.016 against a within-group sd of 0.019 — **the
null still holds**. This comparison is protected because all nine samples sit at the same
count scale, so the transform bias is common to them.

🟡 That is an informative negative. If capture were the GC-biased step, diluting probe 100×
should worsen it as hybridisation moves further from completion. It does not. So on this data
the effect does **not** localise to hybridisation/capture, which pushes it back onto the shared
PCR/sequencing path — consistent with the Schmierer finding that plasmid and genomic DNA give
the same magnitude.

### 🔴 Correction — what the preprint's GC claim actually is

⚠ An earlier version of this file treated the preprint's claim as being about **template
accessibility** (plasmid vs genomic DNA), and then reported that a Schmierer plasmid-vs-genomic
comparison does not support it. That was a misreading.

🟢 The preprint's claim is a comparison of **our two protocols on the same material**:
conventional PCR versus CRISPR-MIP. The accessibility paragraph ("If naked plasmid DNA is more
accessible than gDNA … low GC%, and thus low melting temperature, has an effect in the first
cycles of PCR") is offered as a *proposed explanation* for the protocol difference, not as the
claim itself.

🔴 Two consequences:

1. **The Schmierer plasmid-vs-genomic result does not test the preprint's claim.** It tests the
   proposed mechanism, on someone else's chemistry. It remains a useful negative for that
   mechanism, but it should not be presented as bearing on the PCR-vs-MIP comparison.
2. **Everything analysed here so far is the padlock arm only.** The PCR arm of the same
   experiment was never run. ⚠ So the headline comparison the preprint rests on has not yet
   been reproduced in this repo.

🟡 Also noted: **gDNA was not digested for most runs** — only the GFPg1 trial and
`digested_gDNA_rep1/2`. So the digest is a small side-experiment here, not a standing part of
the protocol, and the "digest as a GC knob" discussion in
[`polymerases.md`](polymerases.md) should be read with that scope.

### The matched PCR-vs-padlock pairs, which is the right test

🟢 From the SDRF, four conditions carry **both** protocols on the same material:

| library | condition | input | PCR runs | padlock runs |
|---|---|---|---|---|
| Brunello-kinome | **`plasmid_lib`** | 0.2 µg plasmid | 4 | 3 |
| Brunello-kinome | `normoxia_3d` | 10.1 µg gDNA | 6 | 4 |
| Brunello-kinome | `normoxia_14d` | 10.1 µg gDNA | 6 | 3 |
| Brunello-kinome | `hypoxia_3d` | 10.1 µg gDNA | 3 | 3 |

🟡 `plasmid_lib` is the cleanest: no cells, no selection, no biology, so any GC difference
between the two protocols is protocol alone.

⚠ The two arms need different parsing, and this is easy to get wrong. The **padlock** guide sits
at a fixed R1 offset 20. The **PCR** guide must be found by **anchor** — the 20 bp immediately
before the scaffold `GTTTTAGAGC` — because the Broad GPP P5 primers carry deliberate **0–8 nt
staggers** for flowcell diversity, so the guide position varies read to read. `extract_pcr.py`
anchors; `extract_mip.py` offsets.

### ✅ Result 2 — the UMI matters, but the first pass of this was wrong

⚠ **Correction.** My first pass took the *median of per-guide logs*. That is unsafe at these
counts: for Poisson, `E[log2 X] = log2(lam) - 1/(2 lam ln2)`, so the transform bias is
**−0.144** at lam = 5 but only **−0.013** at lam = 55. UMIs/guide is ~5–7 and reads/guide ~55,
so the UMI statistic carried ~10× the bias of the read statistic — and because the counts
differ between GC bands, so did the bias. It manufactured part of the trend.

✅ **Redone by aggregating counts within a band first, then taking the log.** Mean ± s.e. over
the nine samples, relative to the 0.45–0.55 band:

| GC band | **reads only** | **UMIs (dedup)** | reads per molecule |
|---|---|---|---|
| 0.00–0.45 | −0.095 ± 0.010 | −0.108 ± 0.007 | +0.013 ± 0.006 |
| 0.45–0.55 | 0 | 0 | 0 |
| 0.55–0.60 | +0.040 ± 0.006 | +0.054 ± 0.005 | −0.014 ± 0.004 |
| 0.60–0.70 | +0.048 ± 0.007 | +0.074 ± 0.005 | −0.026 ± 0.004 |
| **0.70–1.00** | **+0.040 ± 0.007** | **+0.083 ± 0.005** | **−0.044 ± 0.004** |

✅ The decomposition is now **exact** (sums, not medians): `+0.0396 = +0.0833 − 0.0437`.

⚠ What changed: the read-count figure at GC > 0.70 moved from +0.003 to **+0.040** — so the
earlier claim that "read counts look completely flat while molecules say +4.6 %" was an
artefact of the log transform. The honest statement is weaker: **read counts understate
high-GC abundance by ≈0.044 log2 ≈ 3 %, roughly halving the real composition signal rather
than erasing it.** 🟢 The `reads per molecule` column barely moved (−0.042 → −0.044), because
the transform bias largely cancels in a per-guide ratio — which is why that column is the one
to quote.

### ⚠ Sequencing error in the UMI — direction is robust, magnitude is not

🔴 No error-collapse was applied (the authors use `UMIClusterer(threshold=2)`), and with a
13-nt UMI every miscalled read lands on a *fresh* UMI (4¹³ = 67 M, so no collisions). Error
UMIs are therefore **proportional to reads**:

| per-base error | P(UMI miscalled) | implied error UMIs | % of observed |
|---|---|---|---|
| 0.05 % | 0.0065 | 22,375 | 6 % |
| 0.1 % | 0.0129 | 44,615 | **12 %** |
| 0.2 % | 0.0257 | 88,697 | **25 %** |
| 0.5 % | 0.0631 | 217,806 | **61 %** |

✅ But the *direction* of the conclusion survives, and that can be shown rather than hoped.
Writing `U_obs(b) = U_true(b) + p·R(b)` and normalising to the reference band, the observed
UMI trend is a **weighted average of the true molecule trend and the read trend**. Since reads
(+0.040) sit *below* UMIs (+0.083), error contamination drags the estimate *downward*, so the
true molecule trend is **steeper** than measured:

| per-base error | error UMI fraction | true molecules | amplification term |
|---|---|---|---|
| 0 | 0 % | +0.083 | **−0.044** |
| 0.1 % | 12 % | +0.089 | **−0.050** |
| 0.2 % | 25 % | +0.097 | **−0.058** |
| 0.5 % | 61 % | +0.148 | **−0.108** |

🟢 **Every correction makes the amplification bias larger, never smaller.** So −0.044 is a
floor; the sign and the existence of the effect are robust, the magnitude is uncertain by
~2–3×. 🔴 Running the authors' edit-distance-2 collapse is required before quoting a number.

🟡 A second error term pushes the same way. 12.1 % of reads carry a guide absent from the
library and are discarded, and on the Schmierer data that loss was strongly GC-dependent
(Spearman +0.59). High-GC guides therefore lose more reads *and* more UMIs, depressing both
abundance columns — so "+0.083 molecules" is itself a floor. ✅ It largely cancels in the
reads-per-molecule ratio, which is a further reason that column is the trustworthy one.

### ✅ PCR vs CRISPR-MIP, all three matched gDNA conditions, UMI errors collapsed

🟢 The comparison the preprint's claim rests on, now run on every deposited matched pair.
Brunello-kinome, 10.1 µg gDNA, 4M read pairs per sample, UMIs collapsed by the authors'
directional rule at edit distance ≤ 2. log2 at GC > 0.70 relative to the 0.45–0.55 band:

| condition | PCR reads | MIP reads | MIP molecules | **dedup effect** |
|---|---|---|---|---|
| normoxia 3d | −0.135 | −0.130 | −0.151 | **−0.021** |
| normoxia 14d | −0.155 | −0.204 | −0.215 | **−0.010** |
| hypoxia 3d | −0.143 | −0.168 | −0.202 | **−0.034** |
| **mean** | | | | **−0.022 ± 0.010** |

🟢 **The dedup effect replicates: same sign in 3/3 independent conditions, mean −0.022.**
✅ It survives UMI error-collapse, which removed 17–24 % of UMIs — on the `normoxia 3d` pair
the raw effect was −0.027 and the collapsed one −0.021, so **~24 % of it was error artefact
and ~76 % is real.**

⚠ **But it is a UMI effect, not a protocol effect.** PCR reads and MIP reads agree far more
closely with each other than either does with MIP molecules, and the PCR-vs-MIP-molecule gap
(+0.016, +0.060, +0.059) is larger and more variable than the dedup effect itself —
consistent with the two arms differing in ordinary ways (depth, efficiency, sample handling)
rather than in GC response.

> 🟡 **Therefore the defensible claim is:** *CRISPR-MIP's UMI corrects a real, reproducible
> GC-dependent amplification bias of ≈0.02 log2 (≈1.5 %) at high GC, which a read-count
> readout cannot correct because it has no molecule counter.* 🔴 The stronger claim — that
> CRISPR-MIP is **more accurate** than PCR — is **not** supported by this data, because
> nothing here measures truth. See [`../EXPERIMENTS.md`](../EXPERIMENTS.md) for the designs
> that would.

### ⚠ Caveats on these numbers

🔴 Error-collapse still not applied — see the sensitivity table above for what it would do.
🟡 4M read pairs per sample gives only ~5–7 UMIs per guide on a 77,441-guide library, so
per-guide statistics are thin and the band medians carry the weight.
🟡 Read layout was verified three ways before use: the construct predicts guide at R1 offset 20
(EXT_ARM 19 + the Pol III +1 G), an empirical offset scan gives 20 with a 200× margin over the
next-best offset, and 12.1 % of reads fall outside the library, a sensible error rate.

✅ So the strongest single analysis available from existing data is **`t0` vs `normoxia_3d`**
at matched reads-per-UMI: same library, same template type, same protocol, 10× input, three
days apart. And separately, **digested vs non-digested at fixed 7 µg** tests the accessibility
mechanism with input held constant — the two halves of the hypothesis, cleanly separated.

🔴 Neither can be run from public data yet: the raw sequencing is still "will be deposited on
ArrayExpress". Both need the local FASTQ or the per-sample count tables.

🔴 No public dataset has yet been identified that combines a sparse UMI, a no-selection
baseline, and a cycle-number contrast in one experiment. If one is needed, BioGRID ORCS is
the index to search, and the criteria are listed in [`README.md`](README.md).
