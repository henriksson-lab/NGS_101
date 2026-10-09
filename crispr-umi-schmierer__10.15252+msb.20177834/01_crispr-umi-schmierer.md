# CRISPR-UMI (Schmierer) — a UMI cloned into the guide library

> 🟢 verbatim from the paper · 🟡 derived · ✅ computed from the real lentiGuide-Puro map
> and asserted in `tools/selftest.py` · 🔴 not published.

**Schmierer B, Botla SK, Zhang J, Turunen M, Kivioja T, Taipale J.** "CRISPR/Cas9 screening
using unique molecular identifiers." *Mol Syst Biol* 2017;**13**(10):945.
doi:[10.15252/msb.20177834](https://doi.org/10.15252/msb.20177834) · PMC5658704 · CC BY

---

## 1. ⚠ Two different things are called "CRISPR-UMI"

| | **Schmierer (this page)** | **CRISPR-MIP** ([`../crispr-mip__10.1101+2024.03.28.587082/`](../crispr-mip__10.1101+2024.03.28.587082/01_crispr-mip.md)) |
|---|---|---|
| Where the UMI lives | **cloned into the library plasmid** | **on a capture probe**, applied afterwards |
| What it labels | a **transduced cell lineage** — one virus, one RSL | a **captured genomic molecule** |
| Counting it gives | how many **clones** survived | how many **molecules** were there |
| Changes the plasmid? | **yes** | no |
| Read out by | PCR off gDNA; the RSL sits in the i7 position | padlock capture, then PCR off the circle |

Both are "CRISPR + UMI" and both improve a screen's statistics, but they intervene at
opposite ends of the experiment. Mixing up which one a vector belongs to is an easy and
expensive mistake.

Schmierer calls the UMI a **Random Sequence Label (RSL)**, which is the less ambiguous name.

## 2. What it is for 🟢

> "Random Sequence Labels (RSLs) are incorporated into the guide library, which act as unique
> molecular identifiers (UMIs) to allow massively parallel lineage tracing and lineage
> dropout screening. RSLs greatly improve the reproducibility of results by increasing both
> the precision and the accuracy of screens. They reduce the number of cells needed to reach
> a set statistical power."

Each virion carries one guide *and* one random 6-bp RSL, so every transduced cell founds a
lineage with a unique guide+RSL tag. Three analyses then become possible from one screen 🟢:

| | | |
|---|---|---|
| **TCA** | total count analysis | conventional read counting, ignoring the RSL |
| **LDA** | lineage dropout analysis | count *distinct RSLs* per guide — how many clones survived, not how many reads |
| **IRA** | internal replicate analysis | treat RSL sublineages within one guide as replicates generated *inside* a single screen |

🟡 The point of LDA: read count conflates "many surviving clones" with "one clone that
expanded". Counting distinct RSLs separates them, which is why precision improves without
more cells.

## 3. The vector: pLenti-Puro-AU-flip-3xBsmBI 🟢

Built from **lentiGuide-Puro (Addgene #52963)** by

> "replacing the sequence `gttttagagctagaaatagcaagttaaaa……TTTTTT` with
> `gtttAagagctagaaatagcaagttTaaa……TTTTTTcgtctct` to create an AU-flip (Chen *et al*, 2013)
> and an additional BsmBI site downstream of the tracrRNA."

**The AU-flip is exactly two substitutions** ✅ — position 5 `T`→`A` and position 26 `A`→`T`:

```
original   GTTTTAGAGCTAGAAATAGCAAGTTAAAAT...
AU-flip    GTTTAAGAGCTAGAAATAGCAAGTTTAAAT...
               ^                    ^
```

🟡 It removes a `TTTT` stretch that Pol III can read as a terminator, and restores base
pairing in the lower stem — the same motivation as the F+E scaffold, by a smaller edit.

> 🔴 **Consequence worth stating loudly:** `GTTTTAGAGCTAGAAATAGCAAG` — the diagnostic head of
> the original scaffold, and the annealing site of several published primers and probes — is
> **absent from this vector** ✅. Anything designed against the scaffold's 5' half silently
> fails. Primers that anneal to the scaffold's *3'* end (`…GGCACCGAGTCGGTGC`) are unaffected ✅.

## 4. The library insert 🟢

```
ggctttatatatcttgtggaaaggacgaaacaccg [20-nt spacer] gtttaagagctagaaatagcaagtttaaa…ttttttt
GATCGGAAGAGCACACGTCTGAACTCCAGTCAC nnnnnn aagcttggcgtaactagatcttgagacaaa
```

| part | what it is |
|---|---|
| `ggctttatat…acgaaacaccg` | U6 3' end through the Pol III **+1 G** — the Gibson overlap |
| `[20-nt spacer]` | the guide, from the array |
| `gtttaagagct…ttttttt` | **AU-flip scaffold** + terminator |
| `GATCGGAAGAGCACACGTCTGAACTCCAGTCAC` | 🟢 the **Illumina Read-2 / Index-1 adapter, built into the construct** |
| `nnnnnn` | the **6-bp RSL** |
| `aagcttggcg…gagacaaa` | lentiGuide-Puro sequence past the HindIII site |

🟡 Note where the RSL lands: immediately 3' of the Illumina Index-1 adapter — i.e. **exactly
the i7 index position**, so it is read by the Index 1 read rather than needing its own
sequencing primer.

**Assembly** 🟢: the array oligo is annealed to a single 119-bp oligo carrying the RSL and
the Illumina site, double-stranded with outer primers, and cloned by **Gibson assembly**
(NEBuilder HiFi), then electroporated into Endura cells and grown at 30 °C.

**Library** 🟢: 2,325 genes, **23,279 guides**, 101 non-targeting controls; guide sequences
taken from Wang *et al.* 2014; synthesised on CustomArray.

## 5. Why this changes the readout ✅

Because the Illumina adapter is already in the construct, the readout is **not** the
two-step PCR used with lentiCRISPRv2/lentiGuide, where PCR1's reverse primer grafts a
Read-2 site on as a non-templated tail. Here that site is templated.

Rebuilding the vector from the real lentiGuide-Puro map and applying the published edits:

| | plain lentiGuide-Puro | pLenti-Puro-AU-flip-3xBsmBI |
|---|---|---|
| cloned vector | 8,323 bp | **8,362 bp** (+39) |
| PCR1 with `CRISPR_PCR1-F/R` | 557 bp | **596 bp** |
| original scaffold head | present | **absent** |
| Illumina Read-2 adapter | absent | **present** |
| CRISPR-MIP ligation arm | absent | absent |

### The published readout: three nested PCRs to 288 bp 🟢

The Methods *do* give the primers. Verbatim from "Library preparation and sequencing":

| | sequence | cycles |
|---|---|---|
| `PCR1_FW` | `GGACTATCATATGCTTACCGTAACTTGAAAGTATTTCG` | 14 |
| `PCR1_REV` | `CTTTAGTTTGTATGTCTGTTGCTATTATGTCTACTATTCTTTCC` | |
| `PCR2_FW` | `TCTTTCCCTACACGACGCTCTTCCGATCTCTTGTGGAAAGGACGAAACAC` | 19 |
| `PCR2_REV` | `AGAAGACGGCATACGAGATCTGCCATTTGTCTCAAGATCTAGTTAC` | |
| `PCR3_FW` | `AATGATACGGCGACCACCGAGATCTACAC` + **[i5]** + `TCTTTCCCTACACGACGCTCTTCCG` | 14 |
| `PCR3_REV` | `CAAGCAGAAGACGGCATACGAGATCTGCCATTTG` | |
| `CRIPSRSEQ` | `CGATCTCTTGTGGAAAGGACGAAACACCG` | custom read primer |

(`CRIPSRSEQ` is the paper's own spelling.) 200 µg of gDNA is split across **40 parallel PCR1
reactions** of 5 µg each, pooled after 14 cycles; PCR2 runs off 5 µl of the pool, PCR3 off
2 µl of PCR2. KAPA HiFi HotStart throughout. The product is gel purified.

✅ Rebuilding the vector and running all three PCRs in `tools/selftest.py` gives **288 bp** —
exactly the figure the paper quotes. That is the strongest single validation of the
reconstruction, and it is what pinned down the terminator length (see the note below).

**Sequencing** 🟢: `20 + 6 + 6` cycles on a HiSeq4000 —

| read | cycles | what it reads |
|---|---|---|
| Read 1, primed with `CRIPSRSEQ` | 20 | the **guide**, from cycle 1 |
| Index 1 (**i7**) | 6 | the **RSL** |
| Index 2 (**i5**) | 6 | the sample index |

⚠ Note the assignment: **i7 reads the RSL, i5 reads the sample index** — the opposite of the
usual convention, and an easy way to lose the lineage labels in demultiplexing.

✅ `CRIPSRSEQ` ends `…GACGAAACACC` + the Pol III **+1 G**, so its 3' end abuts the spacer with
nothing in between; Read 1 covers the whole 20-nt guide and no leader has to be skipped.
This is what Fig 1 means by "a custom primer placed directly upstream of the guide".

### ✅ The same readout runs on plain lentiGuide-Puro — 39 bp shorter

Running the identical three PCRs on a cloned plain lentiGuide-Puro map gives **249 bp**, not
288. The difference is **39 bp = the 33-nt Illumina adapter + the 6-nt RSL**, i.e. precisely
what Schmierer inserted. So:

> **A 288 bp vs 249 bp PCR3 product distinguishes the two vectors outright** — no sequencing
> needed. Both templates amplify, which is exactly why the mix-up is easy to miss.

🟡 On lentiCRISPRv2 the readout fails at PCR2 instead: `PCR2_REV`'s site is absent, because
the sequence downstream of the terminator differs. Only the lentiGuide lineage carries it.

### ✅ Correction found by this check: the terminator is 6 T, not 7

Reproducing 288 bp required the U6 terminator to be **6 T**, which is what the real
lentiGuide-Puro map carries after the scaffold (`…GTGC` + `TTTTTT` + `AAGCTT…`). An earlier
draft of this page used 7 T, which inflated every downstream figure by 1 nt. All derived
numbers on this page and in `tools/crisprumi.py` have been corrected.

## 6. A padlock probe that works on this vector ✅

Yes, one can be designed — and it turns out *better* than the original, because the capture
can be made to include the RSL as well as the guide. Designed with
`tools/design_padlock.py`, verified end to end against the rebuilt vector.

### The trap to design around 🔴

The CRISPR-MIP probe backbone carries `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC`. This vector
carries `GATCGGAAGAGCACACGTCTGAACTCCAGTCAC` — **a substring of it**. ✅ So a capture that
spans the vector's built-in adapter puts **two copies of the Read-2 primer site into one
amplicon**, which no readout can disambiguate. That constraint, not the AU-flip, is what
actually shapes the design.

### ⚠ The RSL is 141 nt away, and that cannot be shortened

✅ Measured from the published insert, counting the Pol III +1 G as position 1: the RSL
occupies positions **137&ndash;142**. Between the guide and the RSL sit **115 nt of fixed
sequence** — the 76-nt scaffold, the 6-nt terminator and the 33-nt built-in Illumina adapter. The extension
arm must be upstream of the guide (in U6) for the guide to be captured at all, so

> **141 nt is the floor for any single padlock that captures both the guide and the RSL.**
> No arm placement reduces it.

For comparison, the published CRISPR-MIP probe fills **112 nt**. Capturing the RSL is
therefore **1.26&times; the only gap length this chemistry has been shown to work at**.

🟡 That is within the range MIP designs normally quote (~100&ndash;250 nt), so it is not
unreasonable *a priori* — but 🔴 **there is no efficiency data for this probe at 141 nt**, and
padlock gap-fill efficiency and uniformity both degrade as the gap grows. Treat the design
below as a candidate to test, not as a drop-in. The guide-only alternative needs a 66-nt gap,
*shorter* than the validated 112, and carries no such risk.

### Option A: capture the guide *and* the RSL (141-nt gap — test before trusting)

```
extension arm   GTGGAAAGGACGAAACACC          unchanged -- U6 is untouched by the AU-flip
ligation arm    AAGCTTGGCGTAACTAGATCTTGAGA   26 nt, Tm 58.5 C
```

✅ Verified: **one capture site**, gap **142 nt**, circle **271 nt**, and the captured region
contains the `+1 G`, the **20-nt spacer**, the full AU-flip scaffold, the terminator, the
vector's Illumina adapter **and the 6-bp RSL**.

Tm margin is 4.4 °C above the extension arm — the same margin the published probe uses, and
in the right direction.

**Two consequences, both good:**

- ✅ The existing `P7_tracrRNA_rev` primer site still lies inside the capture, so that primer
  is unchanged.
- 🟡 The probe backbone **must drop its own Read-2 site**, letting the vector's serve instead.
  The RSL then sits in the index position and is read as the i7 — which is exactly what
  Schmierer designed it to do. The probe's own UMI still counts molecules, so you get
  *molecule* counts and *lineage* labels from one amplicon.

### Option B: guide only — shorter gap than the validated probe

```
ligation arm    TCAACTTGAAAAAGTGGCACCGA      23 nt, Tm 58.3 C
```

Sits in the scaffold's 3' half at positions 45–68, clear of both AU-flip positions (5 and 26).
✅ Verified to give **exactly one capture, gap 66, circle 192 nt, spacer captured** on
**lentiCRISPR v1, lentiCRISPRv2, lentiGuide-Puro and the Schmierer vector alike**.

Trade-offs: it does **not** capture the RSL, and ✅ the 66-nt capture is too short to hold
both existing PCR primer sites, so that pair needs redesigning.

🟡 **If the 141-nt gap turns out not to work**, this is the fallback — but note what is lost:
the guide and the RSL are then no longer in the same molecule, so a padlock readout gives
molecule counts per guide while the lineage labels have to come from the native PCR readout
separately. The pairing of guide to lineage, which is the whole point of the RSL, is not
recoverable from two independent assays.

### Summary

| | A: guide + RSL | B: guide only |
|---|---|---|
| ligation arm | `AAGCTTGGCGTAACTAGATCTTGAGA` | `TCAACTTGAAAAAGTGGCACCGA` |
| gap / circle | 142 / 271 nt | **66 / 192 nt** |
| vs the validated 112-nt gap | **1.27&times; — untested** | 0.59&times; — comfortably shorter |
| captures the guide | yes | yes |
| captures the RSL | **yes** | no |
| works on other vectors | no | **v1, v2, lentiGuide, this one** |
| existing PCR primers | P7 unchanged; backbone must lose its Read-2 site | both need redesign |
