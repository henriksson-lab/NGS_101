# The barcode cassette — what is known, and a hypothetical model

The single real gap in the construct. Everything here that is not marked 🟢 is **hypothesis**,
and must be drawn on the page as such.

---

## 1. Two different things are both called "adapter"

| Component | P/N | Structure known? |
|---|---|---|
| **Ligation Adapter** — the top/bottom oligo pair in the libprep guide | CRP-LGA1 | 🟢 **Fully known and verified.** 12 bp stem, 20 nt single-stranded 3' arm, **1-nt 3'-T overhang** on the bottom strand, `/5Phos/` on top, `/5AmMC6/` on bottom. See `04_nebnext_illumina.md` §5. |
| **Barcode cassettes** — the A/B/C/D plate used in the four split-pool rounds | CRP-PLT1 | 🔴 **Nothing published.** No sequences, no strand design, no overhangs, no modifications. |

So: we *do* have the ligation adapter's overhang. We have **nothing** for the barcode cassettes —
and those are the ones that matter for the left half of the construct.

## 2. Constraints any model must satisfy

These are hard facts the hypothesis has to fit:

1. 🟢 **Input is dA-tailed.** The barcoding kit's end-prep (20 °C 30 min → 65 °C 30 min) is a
   standard blunt + dA-tail profile. So the **round-A cassette must ligate to a 3'-dA end** — i.e.
   it needs a **1-nt 3'-T overhang**, the same TA chemistry as the ligation adapter.

   > ⚠ **But not every PTA amplicon end can be dA-tailed.** A PTA amplicon's 3' terminus is an
   > α-thio-dideoxynucleotide: no 3'-OH, phosphorothioate linkage, exonuclease-resistant. It cannot
   > be extended, A-tailed, or rescued. 🟡 So barcoding must occur at the ends that *do* carry a
   > normal 3'-OH — from polymerase dissociation before terminator incorporation, from nicks, and
   > from fill-in at recessed ends. Each ~1 kb amplicon has two ends and only needs one. The 5' ends
   > additionally carry a **5'-OH** (the synthetic random primer), so the kinase activity in the
   > end-prep mix is load-bearing. See `06_pta.md` §2.
2. 🟢 **Read-2 layout is `8 + 4 + 8 + 4 + 8 + 4 + 8 + 1`** (Bascet). Four 8 nt barcodes, **three**
   4 nt linkers *between* them, then **one** base, then genomic sequence.
3. 🟢 **Ligation order A → B → C → D**, with D outermost; R2 reads D, C, B, A, then insert.
4. 🟢 **The P7 / Read-2 arm cannot come from the ligation adapter** — that adapter installs only the
   Read-1 / P5 side, at both ends (verified in `04_nebnext_illumina.md` §5). The Read-2 arm must
   therefore be installed on the barcode side: either by the round-D cassette or by a terminal
   adapter applied after round D.
5. 🟢 For the Atrandi i7 primer to prime, the molecule must present
   `GTGACTGGAGTTCAGACGTGTGCTCTT` (the primer's 3' 27 nt) on the barcode-distal end.

> **Constraint 2 is quietly informative.** Three linkers for four barcodes, plus a single base
> before the insert, is exactly what TA ligation at the genomic junction predicts: the "+1" is the
> **T/A junction base**, and each of the three barcode-to-barcode junctions contributes a 4 nt
> cohesive overhang. Bascet's source comment — *"8 barcodes, 3 spacers, and 1 to account for
> ligation"* — says precisely this.
>
> ⚠ Atrandi's own Figure 6 instead draws a **fourth** linker at positions 45–48, between barcode A
> and the insert. That contradicts constraint 2. Figure 6 is illustrative (its example sequence
> matches no real barcode), so we follow Bascet. Flagging it only so the choice is visible.

## 3. Hypothetical model

**HYPOTHESIS — not from any Atrandi document.** Each cassette is a short duplex that (a) joins the
growing molecule through a cohesive end and (b) leaves behind a fresh 4 nt overhang for the next
round. All 24 barcodes within a round share the same linker, so one linker sequence per round.

### Round A — joins the dA-tailed genome (TA ligation)

```
        dA-tailed PTA amplicon                 round-A cassette (hypothetical)

  5'- XXXXXXXXXXXX...XXX A                 5'-p AAAAAAAA LLLL -3'
  3'- XXXXXXXXXXXX...XXX                     3'- aaaaaaaa      -5'
                        ^                       ^^^^^^^^ ^^^^
                   3'-dA overhang          8-nt barcode  4-nt 5' overhang
                                                         left for round B
```

> The cassette's 3'-T overhang (on the strand not drawn above for clarity) pairs with the amplicon's
> 3'-dA. **This junction base is the "+1" in Bascet's trim length.**

### Rounds B, C, D — cohesive 4 nt junctions

```
     growing molecule (ends in a 4-nt 5' overhang)      incoming round-B cassette

  5'-       LLLL AAAAAAAA A XXX...                   5'-p BBBBBBBB LLLL -3'
  3'- ...                                              3'- bbbbbbbb      -5'
            ^^^^                                                    ^^^^
       overhang left by round A                       new overhang for round C
```

After four rounds the top strand reads, 5'→3':

```
[Read-2 / P7 arm] DDDDDDDD LLLL CCCCCCCC LLLL BBBBBBBB LLLL AAAAAAAA T XXXXXX...
                  <--8--> <-4-> <--8--> <-4-> <--8--> <-4-> <--8-->  ^   insert
                                                                     TA junction
```

which is exactly the Bascet layout. ✓

### The D-side (Read-2 / P7) arm — a qualified guess 🟡

The outer end of the round-D cassette is **semi-known**: it must be a TruSeq Read-2 arm, because
Atrandi's i7 primer has to land on it. The open question is *how long*, and the primer sizing tells
us something.

**The design rule, established on the P5 side where both halves are known:**

| | |
|---|---|
| Ligation adapter's single-stranded arm | `GTCGTGTAGGGAAAGAGTGT` — **20 nt** |
| Atrandi P5 primer's 3' annealing portion | `ACACTCTTTCCCTACACGAC` — **20 nt** |
| Relationship | exact reverse complement, **identical length, zero slack** |

Atrandi size the primer to the arm exactly. Applying the same rule to the i7 side:

```
Atrandi i7 primer       CAAGCAGAAGACGGCATACGAGAT AACCTG GTGACTGGAGTTCAGACGTGTGCTCTT
                        +------- P7, 24 ------+ +6 nt+ +----- 27 nt annealing -----+

PREDICTED round-D arm (top strand, 5'->3', immediately 5' of barcode D):

                                               GTGACTGGAGTTCAGACGTGTGCTCTT      27 nt

canonical TruSeq Read 2                        GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT   34 nt
                                                                          ^^^^^^^
                                                        i.e. canonical minus the terminal CCGATCT
```

**Tm corroborates the sizing.** Nearest-neighbour Tm (SantaLucia 1998; 0.5 µM primer, 50 mM Na⁺ —
Q5 buffer differs, so read these as relative, not absolute):

| Primer, 3' annealing portion | Length | GC% | Tm |
|---|---|---|---|
| Atrandi **P5** `ACACTCTTTCCCTACACGAC` | 20 | 50.0 | **54.3 °C** |
| Atrandi **i7** `GTGACTGGAGTTCAGACGTGTGCTCTT` | 27 | 51.9 | 63.0 °C |
| NEB/Illumina i5 (full) | 33 | 51.5 | 66.9 °C |
| NEB/Illumina i7 (effectively annealed, 32) | 32 | 53.1 | 66.3 °C |

- **Atrandi anneal at 54 °C; the P5 primer's Tm is 54.3 °C.** The annealing temperature is set by
  the *limiting* primer — the short P5 arm — not by the i7 arm at 63 °C.
- **NEB anneal at 65 °C; their arms give 66–67 °C.** Same logic, longer arms.
- This is why Atrandi cannot simply tell you to use NEB's primers at NEB's cycling conditions: their
  P5-side arm is 13 nt shorter than canonical, so a full-length Illumina P5 primer would have 13
  unpaired 3' bases and could not prime at all, and the whole reaction has to drop to ~54 °C.

**Confidence, stated honestly.** The P5-side truncation is *proven* — the adapter sequence is
published and the arithmetic is exact. The i7-side arm length is a **qualified guess by analogy**:

- *Supporting it:* the demonstrated primer-sized-to-arm design rule; the fact that Atrandi print
  their own i7 primer rather than only citing NEB's; and the consistent 54 °C/65 °C split.
- *Against it:* a short primer is also compatible with a **full-length 34 nt arm** — the primer
  would simply land 7 nt further out and still extend normally. And Atrandi may print their own
  primer merely because they use a **6 nt single index** with no i5, which is format-incompatible
  with NEB's 8 nt dual-index sets — a simpler explanation that has nothing to do with arm length.

**Therefore: draw 27 nt, label it as inferred, and note the 34 nt alternative in the caption.**
Either way the primer-binding site is correct and the construct is drawable; only the 7 nt
`CCGATCT` between the arm and barcode D is uncertain.

### Where the P7 arm enters — two candidate models

| | Model | Consequence |
|---|---|---|
| **M1** | The **round-D cassette** itself carries the Read-2 / P7 arm on its distal side. | Simplest; no extra step. The distal end must be ligation-dead (blocked or non-phosphorylated) so the ligation adapter cannot attach there during libprep. |
| **M2** | A **separate terminal adapter** is ligated after round D. | Would explain why `Ligation Adapter (CRP-LGA1)` ships in the *barcoding* kit while being used in the *libprep* guide. |

🟡 **M1 is the better fit.** Whichever is true, the barcode-distal end **must be unavailable for
ligation** during libprep — otherwise it would receive the Read-1/P5 arm like every other end, the
molecule would carry P5 at both ends, and nothing would amplify exponentially. This is the same
suppression-PCR logic described in `04_nebnext_illumina.md` §5.

## 4. How to draw this on the page

- Give the cassette its own `<h3>` step with an explicit caption in the `<i>…</i>` style the
  scg_lib_structs pages use for sub-products:

  > *Hypothetical — Atrandi do not publish the barcode cassette design. The 4 nt cohesive overhang
  > shown here is inferred from the read layout (8+4+8+4+8+4+8+1) and from the dA-tailed input; the
  > real chemistry may differ.*

- Use placeholders throughout, per the agreed convention: `AAAAAAAA` / `BBBBBBBB` / `CCCCCCCC` /
  `DDDDDDDD` for barcodes, `LLLL` for linkers.
- Draw the overhang explicitly rather than hiding it — the point of the step is to show *where the
  joint is*, not to assert its sequence.
- Carry a one-line caveat in the `<info>` preamble so a reader never mistakes the cassette for
  documented chemistry.

## 5. What would settle it

- Baronas et al., *Science* (2026), doi:[10.1126/science.ady7227](https://doi.org/10.1126/science.ady7227)
  — the Atrandi SPC platform paper, the most likely public home of the oligo architecture.
- Atrandi directly.
- Reading the linker bases straight off our own R2 data at offsets 8–11, 20–23 and 32–35: if the
  three linkers are constant across reads, that both confirms the model and hands us the sequences.
