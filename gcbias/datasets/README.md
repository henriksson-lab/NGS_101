# Dataset catalogue — public data for validating PCR GC bias

Why this exists: reviewers of **CRISPR-MIP** asked for more on PCR GC bias, and the plan is
to validate it against public data. Suitable data turns out to be scarce, and the reason is
specific and worth writing down once rather than rediscovering.

Marking convention, as elsewhere in this repo:
🟢 verbatim from the paper · 🟡 derived or inferred · 🔴 not published · ✅ computed here.

## The documents

| | paper | what it has |
|---|---|---|
| [`crispr-umi-schmierer.md`](crispr-umi-schmierer.md) | Schmierer 2017, *Mol Syst Biol* 13:945 | 5 ENA runs, RKO essentiality screen, 6-nt RSL |
| [`crispr-umi-michlits.md`](crispr-umi-michlits.md) | Michlits 2017, *Nat Methods* 14:1191 | 4 SRA runs / 9 original BAMs, mESC + MEF screens, 10-nt barcode |
| [`crispr-mip.md`](crispr-mip.md) | our preprint, 2024.03.28.587082 | own data (not yet deposited) + the two external sets already used |
| [`crispr-star.md`](crispr-star.md) | Ueberschär 2025, *Nat Biotechnol* | GEO GSE262309 — 10-nt UMI **in a read**, already as sgRNA-UMI count tables |
| [`polymerases.md`](polymerases.md) | all four | which enzyme, buffer and cycling each paper used at each amplifying step |

## Headline verdict

Both CRISPR-UMI papers have a lineage UMI. **Neither gives you a usable molecule counter in
its public data** — but for two completely different reasons, and they are nearly opposite:

| | Schmierer | Michlits |
|---|---|---|
| UMI length | 6 nt (RSL) | 10 nt |
| UMI in the public FASTQ? | ✅ **yes** — in the index field of the FASTQ header, but only in the *author-submitted* files | ✅ **no** — `nreads=1`, the index reads are not in the SRA run |
| UMI space vs usage | ✅ **81.9 % occupied — saturated** | ✅ 0.300 % occupied — sparse |
| lineages distinguishable | ✅ ~68 % | ✅ ~99.9 % |
| usable for dedup? | the UMI is there but **too degenerate** | the UMI is clean but **not accessible** |

So: the dataset whose UMI you can read is the one whose UMI does not resolve molecules, and
the dataset whose UMI would resolve molecules is the one whose UMI was not submitted.

✅ **[CRISPR-StAR](crispr-star.md) (2025) appears to break that deadlock** — same 10-nt
barcode as Michlits, but carried in read 2 rather than a stripped index, and GEO serves
per-sample sgRNA-UMI count tables directly. It also adds a different sequencing platform and
a built-in active/inactive internal control. The open question there is joining its
`<gene>_<n>` guide names to sequences; see that page.

## ⚠ The UMI only means anything paired with its guide

A bare 6-nt RSL has **4,096** possible values. Against a 23,250-guide library that is
saturated by construction — it carries no molecule information on its own. It works only
because the RSL is read *in the same molecule as its guide*, making the unit of counting the
**(guide, UMI) pair**, not the UMI. All occupancy numbers below are therefore computed over
the paired space, which is the only space that is meaningful.

✅ Even paired, Schmierer is nearly full:

```
                      guides   UMI space   possible pairs    observed      occupancy
Schmierer  6-nt RSL   23,250       4,096      95,232,000    78,000,000 🟢    81.9 %
Michlits  10-nt BC    26,514   1,048,576  27,801,944,064    83,500,000 🟢     0.300 %
```

✅ Schmierer's paired space is **273× more saturated** than Michlits'. Within one guide,
~3,355 lineages are drawn from 4,096 labels, so only about **68 %** of them land on a label
no other lineage in that guide shares. Michlits draws ~3,149 lineages from 1,048,576 labels:
**99.9 %** separable.

> ⚠ **Do not naively deduplicate the Schmierer data.** This is also what the authors say in
> person: the lineage UMI is not very unique. Collapsing by distinct (guide, RSL) pairs there
> does not count molecules — it counts *occupied label slots*, a quantity that saturates
> toward 4,096 per guide and therefore **compresses exactly the abundant guides** you would
> want to compare. Any GC-vs-abundance trend measured that way is partly the saturation
> curve, not chemistry. The paper itself filters low-count RSLs before use 🟢, which changes
> the effective denominator again.

## What would make a dataset suitable

Collected here because the search is the hard part. A dataset is usable for separating
**PCR/amplification bias** from **biology** if it has:

1. **UMI physically linked to the guide**, both readable from the same molecule — otherwise
   there is no (guide, UMI) pair and see the warning above.
2. **UMI space far larger than lineages per guide.** Occupancy well under ~10 %; otherwise
   distinct-UMI counts are a saturating function of abundance and confound the very axis
   being measured.
3. **The UMI actually present in the public files.** Index reads are routinely dropped at
   deposition — this is the single most common disqualifier, and it is invisible in the
   metadata. Check `nreads` and look at a real record.
4. **Known library sequences**, so GC% per guide can be computed.
5. **A low-selection sample** — plasmid input or a very early timepoint — so a GC trend
   cannot be explained by differential growth.
6. **Contrast in amplification**: different cycle numbers, or an enrichment step that lets
   cycles be cut (which is CRISPR-MIP's and Michlits' own argument).
7. **Raw reads, not pre-collapsed count tables** — a processed matrix has already made the
   dedup decision for you, usually undocumented.

🟡 Items 2 and 3 are what rule out the two CRISPR-UMI papers, and they are the reason this
is hard. Items 5 and 6 are what make the *plasmid-input* samples disproportionately valuable
in both: no biology, and in Michlits' case a deliberately low cycle count.

## The analysis built on this

[`..`](../README.md) holds the download scripts, the molecule-count
model, and the GC-bias comparison. It also carries the result that matters for item 2 above:
naive deduplication does not merely lose power, it **manufactures** a GC slope, and the size
of the artefact is calibrated by simulation there.

## Re-deriving the numbers

```
python3 ref/datasets/tools/selftest.py
```
Pins the occupancy arithmetic, the accessions, and the read-structure facts, so a later edit
cannot quietly change a figure quoted in these documents.
