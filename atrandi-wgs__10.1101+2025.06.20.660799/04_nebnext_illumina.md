# NEBNext + Illumina oligo chemistry — sequence reference

Everything needed to draw the library construct base-by-base on the NEBNext/Illumina side.
All sequences below were read out of the vendor's own current PDF; the derived claims
(reverse complements, pairing registers, lengths) were verified computationally.

**Confidence marks:** 🟢 verbatim from the vendor PDF · 🟡 derived/verified arithmetic ·
🔴 not published.

---

## 1. Canonical Illumina sequences 🟢

Source: **Illumina, "Illumina Adapter Sequences", Document # 1000000002694, v22, September 2025** —
<https://support-docs.illumina.com/SHARE/AdapterSequences/1000000002694_22_ilmn-adapter-sequences.pdf>

```
P5                     AATGATACGGCGACCACCGAGATCTACAC            29 nt
P7                     CAAGCAGAAGACGGCATACGAGAT                 24 nt
TruSeq Read 1 primer   ACACTCTTTCCCTACACGACGCTCTTCCGATCT        33 nt
TruSeq Read 2 primer   GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT       34 nt
```

Adapter-trimming (read-through) sequences, both strands:

```
Seen in Read 1 (= the Read-2-side arm)   AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
Seen in Read 2 (= the Read-1-side arm)   AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
```

🟡 Verified: `revcomp(Read 1 primer) = AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` exactly;
`revcomp(Read 2 primer) = AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC`. **The leading `A` in both is the
dA tail, not part of the adapter.** The shared 12 bp duplex core `GATCGGAAGAGC` / `GCTCTTCCGATC` is
the *same* on both arms — so calling `GATCGGAAGAGC` "the Read-2 stem" is convention only.

⚠ Illumina's document never uses the strings "P5" or "P7" — those are community convention for the
flow-cell lawn oligos. What the document prints is `TruSeq Universal Adapter`, `Index 1 (i7)
Adapters`, `Index 2 (i5) Adapters`, etc.

## 2. NEBNext Ultra II FS DNA Library Prep Kit (NEB #E7805S/L, #E6177S/L) 🟢

Manual `manualE6177-E7805.pdf`, **Version 4.0_7/23** —
<https://www.neb.com/en/-/media/nebus/files/manuals/manuale6177-e7805.pdf>

> **The kit contains no oligos at all.** Manual p.1: *"NEBNext Oligo Kit options can be found at
> www.neb.com/oligos. Alternatively, customer supplied adaptor and primers can be used."*
> This is the door Atrandi — and our protocol — walk through.

Kit components: E7806A FS Enzyme Mix · E7807A FS Reaction Buffer · E7648A Ultra II Ligation Master
Mix · E7374A Ligation Enhancer · E7649A Ultra II Q5 Master Mix · E7808A TE.

### Fragmentation + end repair + dA-tailing (one tube, no cleanup)

| Reagent | µL |
|---|---|
| DNA (100 pg – 100 ng) | 26 |
| FS Reaction Buffer | 7 |
| FS Enzyme Mix | 2 |
| **Total** | **35** |

Heated lid 75 °C: **5–30 min @ 37 °C** (fragmentation) → **30 min @ 65 °C** (end repair + dA-tailing
+ inactivation) → 4 °C.

| Target fragment size | 37 °C time | Optimisation range |
|---|---|---|
| 100–250 bp | 30 min | 30–40 min |
| 150–350 bp | 20 min | 20–30 min |
| 200–450 bp | 15 min | 15–20 min |
| 300–700 bp | **10 min** | 5–15 min |
| 500 bp–1 kb | 5 min | 5–10 min |

> Atrandi specify exactly 10 min → **300–700 bp**.

### Ligation + USER

| Reagent | µL |
|---|---|
| FS reaction mixture | 35 |
| Ultra II Ligation Master Mix | 30 |
| Ligation Enhancer | 1 |
| NEBNext Adaptor | 2.5 |
| **Total** | **68.5** |

20 °C, 15 min, **lid off** → add **3 µL USER Enzyme** → 37 °C, 15 min, lid ≥47 °C.

Adaptor working concentration: 15 µM for 100–500 ng input; **1.5 µM (1:10) for 5–99 ng**;
0.6 µM (1:25) for <5 ng.

> USER is only needed for the non-indexed hairpin adaptor. With the indexed UMI adaptor, or with
> Atrandi's single-tailed adapter, it is omitted.

### PCR

| Reagent | µL |
|---|---|
| Adaptor-ligated DNA | 15 |
| Ultra II Q5 Master Mix | 25 |
| i7 primer | 5 |
| i5 / universal primer | 5 |
| **Total** | **50** |

98 °C 30 s ×1 → [98 °C 10 s / **65 °C 75 s**] ×3–13 → 65 °C 5 min → 4 °C. Heated lid 105 °C.

> 🟢 Critical footnote: *"NEBNext adaptors contain a unique **truncated** design. Libraries
> constructed with NEBNext adaptors require a **minimum of 3 amplification cycles** to add the
> complete adaptor sequences for downstream processes."*

## 3. The NEBNext hairpin Adaptor 🟢

**Not in the E7805 manual** — printed only in the oligo-kit manuals (#E7780 v3.0_6/24 p.2; identical
in #E7600 v6.0_6/24 and #E7335 v8.0_6/24). Component **#E7601A**, supplied at 15 µM.

```
5'-/5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTC[dU]ACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3'
```

`-s-`/`*` = phosphorothioate. Total **65 nt**.

| nt | Sequence | Role |
|---|---|---|
| 1–12 | `GATCGGAAGAGC` | stem, 5'-phosphorylated half |
| 13–31 | `ACACGTCTGAACTCCAGTC` | loop, Read-2-side arm — **truncated**, 2 nt short of `…AGTCAC` |
| **32** | **`dU`** | the uracil, exactly mid-loop |
| 33–52 | `ACACTCTTTCCCTACACGAC` | loop, Read-1-side arm |
| 53–64 | `GCTCTTCCGATC` | stem, complementary half |
| 65 | `T` (phosphorothioate) | 3'-T overhang |

🟡 `revcomp(GCTCTTCCGATC) = GATCGGAAGAGC` → stem is exactly **12 bp**, loop **40 nt**.

```
TOP   5'-p G A T C G G A A G A G C --> A C A C G T C T G A A C T C C A G T C -+
           | | | | | | | | | | | |                                          [dU]   40-nt loop
BOT   3'-T C T A G C C T T C T C G <-- C A G C A C A T C C C T T T C T C A C A -+
          ^+----- 12-bp stem -----+
          3'-T overhang
```

**USER** (🟢 E7780 p.2): *"The loop contains a U, which is removed by treatment with USER Enzyme
(a combination of UDG and Endo VIII), to open up the loop and make it available as a substrate for
PCR."* UDG excises the uracil base; Endo VIII cleaves the backbone at the abasic site. 🟡 The two
resulting arms:

```
Read-2 arm (carries the 5'Phos):  5'-/5Phos/GATCGGAAGAGCACACGTCTGAACTCCAGTC-3'   31 nt
Read-1 arm (carries the 3'-T):    5'-ACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3'       33 nt
```

🟡 **The "truncated design":** the Read-1 arm is *exactly* the full 33 nt TruSeq Read 1 primer, but
the Read-2 arm is 2 nt short of canonical. So in cycle 1 the i7 primer anneals over 32 bp with its
5'-most `GT` hanging off unpaired, and the i5 primer has only 13 bp of 3' complementarity until the
i7 extension product exists — hence the 3-cycle minimum.

## 4. NEBNext Multiplex Oligos — index primers 🟢

### Dual Index Primers Set 2, NEB #E7780S — v3.0_6/24
<https://www.neb.com/en/-/media/nebus/files/manuals/manuale7780.pdf>
Contains the adaptor (#E7601A, 15 µM), USER (#E7602A), i5 primers i509–i516, i7 primers i713–i724,
all 10 µM, 96 reactions.

🟡 Structure: i5 = `AATGATACGGCGACCACCGAGATCTACAC`(29) + **8 nt index** +
`ACACTCTTTCCCTACACGACGCTCTTCCGATC*T`(33) = **70 nt**.
i7 = `CAAGCAGAAGACGGCATACGAGAT`(24) + **8 nt index** +
`GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC*T`(34) = **66 nt**.
The i7 index is carried as the **reverse complement** of what the index read reports.

**i5 (#E7781A–E7788A):**

| Primer | Index in oligo | Index read (revcomp) |
|---|---|---|
| i509 | TTGCTTGC | GCAAGCAA |
| i510 | GAGAGGTT | AACCTCTC |
| i511 | ACCTGGTT | AACCAGGT |
| i512 | AAGCGGAA | TTCCGCTT |
| i513 | CGGAACAA | TTGTTCCG |
| i514 | GGTAAGCT | AGCTTACC |
| i515 | TGTGGCAT | ATGCCACA |
| i516 | ACTACGGA | TCCGTAGT |

**i7 (#E7789A–E7800A):**

| Primer | Index in oligo | Index read (revcomp) |
|---|---|---|
| i713 | AGGAGGAA | TTCCTCCT |
| i714 | AGCAAGCA | TGCTTGCT |
| i715 | TCATCACC | GGTGATGA |
| i716 | CGTAGGTT | AACCTACG |
| i717 | TCAGATCC | GGATCTGA |
| i718 | CGTGATCA | TGATCACG |
| i719 | AGTCGCTT | AAGCGACT |
| i720 | GAACGCTT | AAGCGTTC |
| i721 | TACGCCTT | AAGGCGTA |
| i722 | CTCATCAG | CTGATGAG |
| i723 | TCTTCTGC | GCAGAAGA |
| i724 | GCTGGATT | AATCCAGC |

### Dual Index Primers Set 1, NEB #E7600S — v6.0_6/24
Adaptor is byte-identical to Set 2's. Indices i501–i508 / i701–i712.

| i5 | Index in oligo | i7 | Index in oligo |
|---|---|---|---|
| i501 | TATAGCCT | i701 | CGAGTAAT |
| i502 | ATAGAGGC | i702 | TCTCCGGA |
| i503 | CCTATCCT | i703 | AATGAGCG |
| i504 | GGCTCTGA | i704 | GGAATCTC |
| i505 | AGGCGAAG | i705 | TTCTGAAT |
| i506 | TAATCTTA | i706 | ACGAATTC |
| i507 | CAGGACGT | i707 | AGCTTCAG |
| i508 | GTACTGAC | i708 | GCGCATTA |
| | | i709 | CATAGCCG |
| | | i710 | TTCGCGGA |
| | | i711 | GCGCGAGA |
| | | i712 | CTATCGCT |

### 6-nt single-index kits, NEB #E7335 / #E7500 / #E7710 / #E7730 — v8.0_6/24
Same hairpin adaptor; the i5 side is an **index-less universal primer**:

```
NEBNext Universal PCR Primer for Illumina
5'-AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC*T-3'        57 nt

Index primer format
5'-CAAGCAGAAGACGGCATACGAGAT[6-nt index]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATC*T-3'
```

Set 1 indices 1–12 (oligo → index read): CGTGAT→ATCACG, ACATCG→CGATGT, GCCTAA→TTAGGC,
TGGTCA→TGACCA, CACTGT→ACAGTG, ATTGGC→GCCAAT, GATCTG→CAGATC, TCAAGT→ACTTGA, CTGATC→GATCAG,
AAGCTA→TAGCTT, GTAGCC→GGCTAC, TACAAG→CTTGTA.

> Set 2 (#E7500, indices 13–27) carries **8** bases in the index region but reads only **6** —
> the known TruSeq quirk ("the two bases following the index adapter sequence can vary").

### UMI adaptors — 🔴 **sequence withheld**

DNA Sets 1–4 = **#E7395 / #E7874 / #E7876 / #E7878** (#E7416 is RNA; #E6440 is non-UMI).
Manual v4.0_6/26.

> 🟢 Revision history p.25: *"3.0 … Deleted adaptor sequence in NEBNext Adaptor for Illumina
> Overview because it was inaccurate. The correct adaptor sequence is proprietary. 9/25"*

**NEB published a UMI adaptor sequence in manual versions 1.0–2.0, has since declared it wrong,
removed it, and will not publish the correct one. Do not copy it from an old manual or a
third-party page.** What is stated: a true forked/Y duplex with a single-base T overhang (not a
hairpin, no USER); 8 nt indices pre-installed; 96 unique dual pairs. UMI length is given
inconsistently as 11 nt (v2.0) vs 12 nt (v4.0) — the sequencing table (i7 index read = 8 cycles
without UMI, 20 with) settles it at **12**. Amplification uses a single undisclosed
"NEBNext Primer Mix" — 🔴 no i5/i7 sequences published.

## 5. The Atrandi ligation adapter, resolved 🟡

Sequences 🟢 from `DGPM02323206001` p.4; structure verified base-by-base.

```
Top adapter      /5Phos/GATCGGAAGAGCGTCGTGTAGGGAAAGAGTG*T     32 nt
Bottom adapter   /5AmMC6/GCTCTTCCGATCT                        13 nt
```

### ⚠ It is NOT a forked / Y adapter

The bottom oligo is 13 nt and **all 12** of its 5'-proximal bases pair with the top strand's first
12, leaving **one** single-stranded arm, not two. It is a **12 bp stem with a single 20 nt
single-stranded 3' extension**.

```
              <------ 12-bp stem ------>  <--------- 20-nt single-stranded 3' arm --------->
TOP  5'-/5Phos/ G A T C G G A A G A G C    G T C G T G T A G G G A A A G A G T G T -3'
                | | | | | | | | | | | |
BOT        3'-T C T A G C C T T C T C G /5AmMC6/-5'
             ^
             3'-T overhang = the ligation end
```

Verified:

| Claim | Result |
|---|---|
| `revcomp(GCTCTTCCGATC) == GATCGGAAGAGC` | ✓ 12/12 Watson-Crick |
| `revcomp(GTCGTGTAGGGAAAGAGTGT)` | `ACACTCTTTCCCTACACGAC` = **first 20 nt of TruSeq Read 1 primer** |
| `GCTCTTCCGATCT` | = **last 13 nt of TruSeq Read 1 primer** |
| `"A" + top adapter` | `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT` = **Illumina's Read 2 trim sequence, exactly** |

### What each piece is

| Piece | Identity |
|---|---|
| BOT `GCTCTTCCGATCT` | last 13 nt of the TruSeq Read 1 primer |
| TOP `GATCGGAAGAGC` | its complement — here **not** a Read-2 stem |
| TOP `GTCGTGTAGGGAAAGAGTGT` | revcomp = first 20 nt of the Read 1 primer → the **i5/P5 primer landing site** |

**→ The Atrandi adapter installs only the TruSeq Read 1 / P5 side, at both ends of every fragment.
It carries no Read 2 / P7 arm at all. The P7 / Read 2 side must therefore be installed by the
split-pool barcode cassette.** This is the structural key to the whole construct.

### Ligation geometry

```
      ADAPTER                           dA-TAILED INSERT

 BOT  3'-T -----+                 +----- 5'-p N N N N ...     (insert strand 1)
                |   T:A pair      |
 TOP  5'-p G A T+                 +-- 3'-A N N N N ...        (insert strand 2, dA tail)
          ^                              ^
    ligase seals TOP-5'p to insert-3'A   ligase seals BOT-3'OH to insert-5'p
```

### Why `/5AmMC6/`

🟡 Inference (Atrandi do not explain it). The 5'-amino-C6 sits at the fork point and makes that 5'
end permanently ligation-incompetent. Unlike a bare 5'-OH it **cannot be rescued by the
polynucleotide-kinase activity carried over from the FS end-repair mix** — which is added straight
into the ligation with no cleanup — so adapter-adapter dimers and concatemers cannot form. It also
blocks 5'→3' exonuclease attack on the short oligo.

### The Atrandi indexing primers

```
Index i1 P7:
  CAAGCAGAAGACGGCATACGAGAT  AACCTG  GTGACTGGAGTTCAGACGTGTGCTCTT          57 nt
  +------- P7, 24 nt ------+ +6 nt-+ +--- Read-2 primer, 5' 27 nt -----+

  Illumina canonical:
  CAAGCAGAAGACGGCATACGAGAT[6 bases]GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT    64 nt
                                                              ^^^^^^^  Atrandi drop these 7

P5:
  AATGATACGGCGACCACCGAGATCTACAC  TCTTTCCCTACACGAC                        45 nt
  +--------- P5, 29 nt --------+  +----- 16 nt ---+

  Illumina canonical:
  AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT             58 nt
                                               ^^^^^^^^^^^^^  Atrandi drop these 13
```

- **Index = 6 nt (`AACCTG`), i7 only.** 🟡 `revcomp = CAGGTT` is what appears in the index read on a
  forward-strand instrument.
- **No i5 index** — the P5 primer is the plain universal primer. Atrandi's Figure 2 marks i5 as
  `*Optional`.
- `*` = 3'-terminal phosphorothioate, protecting against Q5's 3' proofreading exonuclease.
- 🟡 **The truncations are not optional.** The P5 primer's 3' 20 nt `ACACTCTTTCCCTACACGAC` are
  *exactly* the revcomp of the adapter's 20 nt arm — a full-length Illumina P5 primer would have
  13 unpaired 3' bases and could not prime in cycle 1. Likewise the i7 primer's 3' 27 nt match the
  first 27 nt of the Read 2 primer, which is what the barcode cassette must present. Same
  "truncated adaptor" trick as NEB, one notch further.
- Consistent with the shorter primers, Atrandi run a **3-step PCR at 54 °C**, not NEB's 2-step 65 °C.

### 🟡 Why the design is a suppression PCR

Because the adapter is symmetric and installs the *same* P5-side arm at both fragment ends, a
fragment carrying adapters but **no** barcode cassette has a P5 landing site at both ends. After one
P5 extension its new 3' end is `AGATCGGAAGAGC`, for which no primer exists — so it amplifies
**linearly**. Only molecules carrying the barcode cassette (and hence the i7 landing site) amplify
exponentially. Almost certainly deliberate.

## 6. Atrandi's modifications to the NEB protocol

| Step | NEB E7805 | Atrandi |
|---|---|---|
| Fragmentation | 5–30 min @ 37 °C, then 30 min @ 65 °C | **exactly 10 min** @ 37 °C ("CRITICAL") |
| Adapter | NEBNext hairpin, 15 / 1.5 / 0.6 µM | **own single-tailed adapter, 2.5 µL of 1.5 µM** |
| Ligase | Ultra II Ligation MM 30 µL + Enhancer 1 µL | identical |
| Ligation | 20 °C, 15 min | identical |
| **USER** | 3 µL, 37 °C 15 min | **omitted — no hairpin to open** |
| Cleanup | 0.8× SPRI | 0.8× AMPure XP |
| PCR enzyme | Ultra II Q5 Master Mix | identical |
| PCR primers | NEBNext i5/i7, 10 µM | 10 µL of 5 µM each; cites #E7780S "or equivalent" but prints its own truncated pair |
| PCR program | 2-step, 65 °C, 3–13 cycles | **3-step, 54 °C anneal, 8–14 cycles** |

## 7. What can and cannot be drawn

| Element | Status |
|---|---|
| NEBNext hairpin adaptor (65 nt, dU position, 5'P, 3'-T, `*`) | 🟢 fully specified |
| NEBNext i5/i7 primers, Sets 1 & 2 (40 sequences) | 🟢 fully specified |
| NEBNext 6-nt single-index primers + universal primer | 🟢 fully specified |
| P5, P7, Read 1, Read 2, both trim strings | 🟢 fully specified |
| Atrandi ligation adapter: sequences, register, geometry | 🟢 sequences, 🟡 register (verified) |
| Atrandi i7/P5 primers, 6-nt index `AACCTG` | 🟢 |
| **NEBNext UMI adaptor** | 🔴 **proprietary; NEB retracted the published sequence** |
| NEBNext UMI kit i5/i7 primers | 🔴 undisclosed "Primer Mix" |
| **Atrandi barcode A/B/C/D oligos + the Read-2/P7 cassette they install** | 🔴 **not published anywhere** |

## References

**NEB**
1. Ultra II FS DNA Library Prep Kit, #E7805S/L, #E6177S/L — Manual **v4.0_7/23** — <https://www.neb.com/en/-/media/nebus/files/manuals/manuale6177-e7805.pdf>
2. Multiplex Oligos, Dual Index Primers Set 2, #E7780S — Manual **v3.0_6/24** — <https://www.neb.com/en/-/media/nebus/files/manuals/manuale7780.pdf>
3. Multiplex Oligos, Dual Index Primers Set 1, #E7600S — Manual **v6.0_6/24** — <https://www.neb.com/en/-/media/nebus/files/manuals/manuale7600.pdf>
4. Multiplex Oligos, Index Primers Sets 1–4, #E7335/#E7500/#E7710/#E7730 — Manual **v8.0_6/24** — <https://www.neb.com/en/-/media/nebus/files/manuals/manuale7335.pdf>
5. Unique Dual Index UMI Adaptors DNA Sets 1–4, #E7395/#E7874/#E7876/#E7878 — Manual **v4.0_6/26** — <https://www.neb.com/en-us/-/media/nebus/files/manuals/manuale7395_e7874_e7876_e7878.pdf>
6. UMI Adaptors DNA Set 1, #E7395 — superseded Manual **v2.0_2/21** (cited only to document the retraction) — <https://www.neb.com/-/media/nebus/files/manuals/manuale7395.pdf>
7. NEB index-combination tool — <https://indexoligo.neb.com/>

**Illumina**
8. **Illumina Adapter Sequences, Document # 1000000002694 v22, September 2025** — <https://support-docs.illumina.com/SHARE/AdapterSequences/1000000002694_22_ilmn-adapter-sequences.pdf>
9. Adapter Sequences landing page — <https://support-docs.illumina.com/SHARE/AdapterSequences/Content/TruSeq-Sequences.htm>
10. *Indexed Sequencing Overview Guide*, Document # 15057455 — forward-strand vs revcomp i5 workflows
11. Illumina Knowledge Article #1800, "Guidelines for reverse complementing i5 sequences for demultiplexing"

**Atrandi** — see `01_barcoding_kit.md` and `02_library_prep.md`.
