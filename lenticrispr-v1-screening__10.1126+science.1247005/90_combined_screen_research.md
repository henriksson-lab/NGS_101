# lentiCRISPR / GeCKO — vectors, cloning and readout

> **Evidence marking.** 🟢 verbatim from a primary source · 🟡 derived · ✅ **computed by us
> from the real Addgene map** and asserted in `tools/selftest.py` · 🔴 not published.

Plasmid maps live in `ref/plasmids/` (see `MANIFEST.md`). Everything marked ✅ is recomputed
on every test run, so it cannot drift.

---

## 1. The vectors

| Vector | Addgene | bp ✅ | GPP name | Role |
|---|---|---|---|---|
| lentiCRISPR v1 | 49535 *(discontinued)* | 13,450 | pXPR_001 | all-in-one Cas9 + guide |
| lentiCRISPRv2 | 52961 | 14,873 | pXPR_023 | all-in-one, ~10× titre |
| lentiGuide-Puro | 52963 | 10,183 | pXPR_003 | guide only |
| lentiCas9-Blast | 52962 | 12,859 | — | Cas9 only |

### ⚠ Two corrections to the common story

**There is no v1-vs-v2 scaffold difference.** ✅ All three guide vectors carry the *identical*
original 76 nt Cong/Ran scaffold. The "F+E" scaffold (`GTTTAAGAGCTATGCTGGAAACAGCATAGCAAG…`)
is absent from all of them; it appears in the Broad GPP vectors (pXPR_011 ✅) and the
Weissman CRISPRi series.

**The real v1→v2 difference that breaks primers is the position of the cPPT/CTS.**

```
lentiCRISPR v1    5'LTR - Psi - RRE - U6-guide-scaffold - [cPPT] - EFS - Cas9-P2A-Puro - WPRE
lentiGuide-Puro   5'LTR - Psi - RRE - U6-guide-scaffold - [cPPT] - EF-1a - Puro        - WPRE
lentiCRISPRv2     5'LTR - Psi - RRE - [cPPT] - U6-guide-scaffold - EFS - Cas9-P2A-Puro - WPRE
                                      ^^^^^^ moved UPSTREAM of U6
```

Every classical reverse readout primer anneals in that cPPT. ✅ On lentiCRISPRv2 the sites
end up in the wrong order on the circle, so the "product" runs the long way round — 12.8 kb,
i.e. nothing usable.

## 2. The BsmBI cloning site ✅

BsmBI = `CGTCTC(1/5)`, 4 nt 5' overhangs. Identical in all three guide vectors:

| | computed |
|---|---|
| BsmBI sites per vector | **2** |
| upstream overhang | **`CACC`** |
| downstream overhang | **`GTTT`** (= the scaffold's own first 4 nt) |
| filler excised | **1,885 bp** |

`revcomp(GTTT)` = `AAAC`, which is exactly the published reverse-oligo overhang. So the
canonical guide oligos fall straight out of the map:

```
top     5'-CACC [G] NNNNNNNNNNNNNNNNNNNN -3'
bottom  5'-AAAC     nnnnnnnnnnnnnnnnnnnn [C] -3'
```

✅ Ligation reconstitutes `…GGACGAAACACC` + **G** + spacer + `GTTTTAGAGCTAGAAATAGCAAG…`
contiguously — verified for both G-initiated and non-G spacers.

> **The +1 G is an off-by-one trap.** `GACGAAACACCG` is the usual search motif, but that
> final G is the Pol III **+1 base**, not part of the promoter. A spacer already starting
> with G supplies it; otherwise it is appended. Treating the motif as "promoter" and then
> adding G again double-counts it — which is why `lib/crispr.py` keeps `U6_3PRIME`
> (11 nt, no G) and `U6_PLUS1` separate. 🟡 GeCKO library spacers are **not** G-initiated
> (first-base composition is roughly uniform), so real library amplicons are 1 nt longer
> than the published figures, which assume a G-initiated spacer. ✅

## 3. Readout primers — the part that differs between protocols

Four published sets. **They are not interchangeable.**

| | **Shalem 2014** | **Zhang F01–F12** | **Joung 2017** | **Broad GPP** |
|---|---|---|---|---|
| Strategy | 2-step nested | 1-step (+adaptor for v2) | **1-step** | **1-step** |
| Fwd anneals | U6 (24 nt) | U6 (24 nt) | U6 (33 nt, no +1 G) | U6 (22 nt) |
| Stagger | 1–9 nt | 1–9 nt | **9–18 nt** ✅ | **0,1,2,3,4,6,7,8** ✅ (no 5) |
| Index position | **inline in read 1** | inline + i7 | i7 | i7 |
| Rev anneals | cPPT | cPPT | **scaffold** | cPPT (KERMIT) *or* EFS (BEAKER) |
| Vector-independent? | no | no | **yes** | no — must pick the right P7 |

### Computed amplicon sizes ✅

On cloned vectors carrying a **G-initiated** 20 nt spacer (add 1 nt for a non-G spacer):

| primer set | lentiCRISPR v1 | lentiCRISPRv2 | lentiGuide-Puro |
|---|---|---|---|
| Shalem F1/R1 (step 1) | 308 bp | **12,760 bp → no product** | 316 bp |
| GPP ARGON + KERMIT | 342 bp | **12,794 bp → no product** | 350 bp |
| GPP ARGON + BEAKER | no product | **285 bp** 🟢✅ | 554 bp → unusable |
| Joung Fwd-1 + KO-Rev-1 | **259 bp** | **259 bp** | **259 bp** |

🟢 GPP publishes "285 nt" for lentiCRISPRv2 — **reproduced exactly**. 🟢 Joung publishes
"~260–270 bp" — we compute 259 (G-spacer) / 260 (non-G). 🟢 GPP publishes "354 nt" for
lentiGuide-Puro; stagger 0 gives 350, so the published figure corresponds to a mid-ladder
stagger.

### ⚠ Picking the right P7 is about product length, not site presence

✅ BEAKER's EFS site **is** present in lentiGuide-Puro (it has its own EF-1α cassette) — it
simply gives a 554 bp product instead of the intended short one. The rule is positional, so
"does my vector contain this site" is the wrong question; "what size product does this pair
give on my vector" is the right one, and `tools/selftest.py` answers it by computing.

| P7 | anneals in | vectors |
|---|---|---|
| **KERMIT** | HIV-1 cPPT/CTS | pXPR_001 (= lentiCRISPR v1), pXPR_003 (= lentiGuide-Puro), and most guide-only vectors |
| **BEAKER** | EFS / EF-1α leader | pXPR_023 (= lentiCRISPRv2) and other all-in-one vectors |
| **MISS PIGGY** | scaffold terminal stem | vector-agnostic |

🟡 "All-in-one ⇒ BEAKER" is **not** a safe rule — several all-in-one pXPR vectors take
KERMIT because, like v1, their cPPT is downstream.

### The v2 adaptor
🟢 lentiCRISPRv2 lacks a cPPT downstream of the guide, so the Zhang lab's workaround grafts
one on: `V2ADAPTOR_R` carries the 25 nt cPPT site as a **non-templated 5' tail** ✅ over a
20 nt EFS-annealing 3' end. One PCR adds the missing site; the standard F01–F12 / R01–R12
primers then work as on lentiGuide-Puro.

## 4. Why every protocol has a stagger

The bases between the sequencing primer and the variable spacer are identical in every
cluster — a monotemplate the sequencer cannot base-call. Five published solutions:

1. **Stagger ladder** — several primers with inserts of different length (GPP, Zhang, Joung, SAM).
2. **Custom Read-1 primer** abutting the cloning junction (Weissman, Bassik, Sabatini — the
   last also needs a custom *index* primer).
3. **Dark cycles** — don't image the invariant bases (CHyMErA: 29 + 20).
4. **Degenerate stagger** in one oligo (GPP P5_MAGNESIUM).
5. **Move the variable region into Read 2 / an mRNA** (all droplet methods, PE/BE sensors).

🟡 The choice is what makes protocols incompatible: ladders are not interchangeable between
labs, and a custom Read-1 primer is **scaffold-specific**, so it cannot be shared across
vector families.

## 5. Library oligo pool 🟢

GeCKO synthesised oligo, 104 nt (from the Zhang lab's own `design_library.py`):

```
TTTCTTGGCTTTATATATCTTGTGGAAAGGACGAAACACCG [20-nt spacer] GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGT
<-- 41 nt: U6 3' end + the +1 G -->                       <-- 43 nt: scaffold nt 1-43 -->
```

Amplified with `Oligo-Fwd` / `Oligo-Knockout-Rev` to a 141 bp insert with **60 bp of
homology on each arm** for Gibson assembly into the BsmBI-cut backbone.

Libraries: GeCKO v2 human **123,411** sgRNAs (A 65,383 + B 58,028; the Addgene PDF's
"122,411" is a typo 🟡), 19,050 genes, 6/gene. Brunello 76,441 sgRNAs, 19,114 genes.

## 6. Sequencing 🟢

Joung: **80 cycles read 1 + 8 cycles index 1**, 5% PhiX (MiSeq) / 20% (NextSeq), >100
reads/sgRNA for plasmid QC and **>500 reads/sgRNA** for a screen. **No custom sequencing
primer** for any of these four sets — all embed the standard TruSeq Read-1 and Index-1
sites, and the stagger supplies the diversity.

Counting: Zhang `count_spacers.py` finds the key `CGAAACACCG` in read positions 30–55 and
takes the next 20 nt. 🟡 That window is tuned for the 9–18 nt Joung staggers; with the
F01–F12 primers the key lands at 23–31, so `KEY_REGION_START` must be set to 0. Broad uses
**PoolQ** with a "search prefix" (`CACCG` for sgRNA) to locate the barcode regardless of
stagger.

## 7. Not published 🔴

GeCKO v1's pooled-library Addgene number; the exact 74 nt v1 array-oligo flank split; stated
amplicon sizes for Shalem's two steps and the F01–F12 set; Shalem's read length; GPP's
read-1 cycle count; cycling conditions for the F01–F12 set.
